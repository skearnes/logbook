# Copyright 2026 Steven Kearnes
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Ligand-only students that try to reproduce the teacher's bioactivity similarity.

Three rungs, cheapest first. Each has to beat the one below it to justify its
extra machinery:

1. ``neighbour`` -- no training at all. A molecule inherits a bioactivity profile
   from its nearest structural neighbours among compounds that *have* been
   assayed. This is the Bioturbo idea (Wassermann et al., JCIM 2013, 53, 692):
   use chemical similarity to map an unannotated molecule into bioactivity
   space, then work there. It is the control the trained models must beat.

2. ``multitask`` -- the standard approach. A network predicts activity across
   all targets at once and its penultimate layer is taken as the fingerprint.
   The profile is supervision, not output.

3. ``distill`` -- trains the actual objective. Rather than predicting the
   profile and hoping similarity follows, it optimises embedding similarity to
   match teacher similarity directly on sampled pairs.

Supervision is masked throughout: ChEMBL records only what was measured, so an
absent compound-target entry means "not tested", never "inactive". Treating
absence as a negative would teach the model that almost everything is inactive
against almost everything, which is both wrong and trivially satisfiable.
"""

from __future__ import annotations

import numpy as np
import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset

import fingerprints

EMBEDDING_DIM = 512
HIDDEN_DIM = 2048
BATCH_SIZE = 256
LEARNING_RATE = 1e-3
WEIGHT_DECAY = 1e-6

DEVICE = torch.device(
    "mps" if torch.backends.mps.is_available()
    else "cuda" if torch.cuda.is_available()
    else "cpu"
)


class Encoder(nn.Module):
    """Structure -> embedding. The embedding is the fingerprint we are after."""

    def __init__(self, n_input: int, n_targets: int):
        super().__init__()
        self.trunk = nn.Sequential(
            nn.Linear(n_input, HIDDEN_DIM),
            nn.ReLU(),
            nn.BatchNorm1d(HIDDEN_DIM),
            nn.Dropout(0.25),
            nn.Linear(HIDDEN_DIM, HIDDEN_DIM // 2),
            nn.ReLU(),
            nn.BatchNorm1d(HIDDEN_DIM // 2),
            nn.Dropout(0.25),
            nn.Linear(HIDDEN_DIM // 2, EMBEDDING_DIM),
        )
        self.head = nn.Linear(EMBEDDING_DIM, n_targets)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.trunk(x)

    def predict_profile(self, x: torch.Tensor) -> torch.Tensor:
        return self.head(torch.relu(self.trunk(x)))


def embed(model: Encoder, features: np.ndarray, batch: int = 4096) -> np.ndarray:
    """Run the encoder over a feature matrix and L2-normalise the output."""
    model.eval()
    chunks = []
    with torch.no_grad():
        for start in range(0, len(features), batch):
            block = torch.from_numpy(features[start:start + batch]).float().to(DEVICE)
            chunks.append(model(block).cpu().numpy())
    out = np.vstack(chunks)
    norms = np.linalg.norm(out, axis=1, keepdims=True)
    return out / np.maximum(norms, 1e-8)


def train_multitask(
    features: np.ndarray, profiles, n_epochs: int = 20, verbose: bool = True
) -> Encoder:
    """Multi-task activity prediction; the penultimate layer becomes the fingerprint.

    ``profiles`` is a sparse matrix of graded activity calls in [-1, 1], where a
    structural zero means "not tested" and is masked out of the loss.
    """
    n_samples, n_targets = profiles.shape
    model = Encoder(features.shape[1], n_targets).to(DEVICE)
    optimiser = torch.optim.AdamW(
        model.parameters(), lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY
    )

    profiles = profiles.tocsr()
    inputs = torch.from_numpy(features).float()

    for epoch in range(n_epochs):
        model.train()
        order = np.random.permutation(n_samples)
        total, seen = 0.0, 0
        for start in range(0, n_samples - BATCH_SIZE + 1, BATCH_SIZE):
            index = order[start:start + BATCH_SIZE]
            # Densified one batch at a time: the full matrix is ~600k x 7k, which
            # would need tens of gigabytes as a dense array.
            y = torch.from_numpy(profiles[index].toarray()).float().to(DEVICE)
            x = inputs[index].to(DEVICE)
            m = (y != 0).float()

            prediction = torch.tanh(model.predict_profile(x))
            # Mean over observed entries only; an unobserved entry contributes
            # nothing rather than being pulled toward zero.
            loss = (((prediction - y) ** 2) * m).sum() / m.sum().clamp(min=1.0)
            optimiser.zero_grad()
            loss.backward()
            optimiser.step()
            total += float(loss) * len(index)
            seen += len(index)
        if verbose:
            print(f"  epoch {epoch + 1:2d}/{n_epochs}  masked MSE {total / seen:.4f}", flush=True)
    return model


def train_distill(
    features: np.ndarray,
    pairs: np.ndarray,
    similarity: np.ndarray,
    n_epochs: int = 20,
    verbose: bool = True,
) -> Encoder:
    """Train embedding cosine similarity to match the teacher's similarity."""
    model = Encoder(features.shape[1], 1).to(DEVICE)
    optimiser = torch.optim.AdamW(
        model.parameters(), lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY
    )

    inputs = torch.from_numpy(features).float()
    pair_tensor = torch.from_numpy(pairs).long()
    target_tensor = torch.from_numpy(similarity).float()

    loader = DataLoader(
        TensorDataset(pair_tensor, target_tensor),
        batch_size=BATCH_SIZE, shuffle=True, drop_last=True,
    )

    for epoch in range(n_epochs):
        model.train()
        total, seen = 0.0, 0
        for pair_batch, target_batch in loader:
            left = inputs[pair_batch[:, 0]].to(DEVICE)
            right = inputs[pair_batch[:, 1]].to(DEVICE)
            target_batch = target_batch.to(DEVICE)

            embedded = model(torch.cat([left, right], dim=0))
            half = len(left)
            predicted = torch.cosine_similarity(embedded[:half], embedded[half:], dim=1)
            loss = nn.functional.mse_loss(predicted, target_batch)

            optimiser.zero_grad()
            loss.backward()
            optimiser.step()
            total += float(loss) * half
            seen += half
        if verbose:
            print(f"  epoch {epoch + 1:2d}/{n_epochs}  pair MSE {total / seen:.4f}")
    return model


def neighbour_profiles(
    query_fps: np.ndarray,
    reference_fps: np.ndarray,
    reference_profiles,
    n_neighbours: int = 25,
    block: int = 256,
) -> np.ndarray:
    """Bioturbo control: inherit a profile from nearest structural neighbours.

    Returns a dense predicted profile per query, the similarity-weighted mean of
    its neighbours' measured profiles. No training involved.
    """
    reference_profiles = reference_profiles.tocsr()
    out = np.zeros((len(query_fps), reference_profiles.shape[1]), dtype=np.float32)

    # The reference side is ~600k x 7k; densifying it whole would need tens of
    # gigabytes, so only the selected neighbour rows are expanded.
    reference = np.ascontiguousarray(reference_fps, dtype=np.float32)
    reference_counts = reference.sum(axis=1)[None, :]

    for start in range(0, len(query_fps), block):
        query = np.ascontiguousarray(query_fps[start:start + block], dtype=np.float32)
        intersection = query @ reference.T
        union = query.sum(axis=1)[:, None] + reference_counts - intersection
        similarity = np.divide(
            intersection, union, out=np.zeros_like(intersection), where=union > 0
        )
        top = np.argpartition(-similarity, n_neighbours, axis=1)[:, :n_neighbours]
        for i in range(len(query)):
            idx = top[i]
            weights = similarity[i, idx]
            total = weights.sum()
            if total <= 0:
                continue
            neighbours = reference_profiles[idx].toarray()
            out[start + i] = (weights[:, None] * neighbours).sum(axis=0) / total
    return out
