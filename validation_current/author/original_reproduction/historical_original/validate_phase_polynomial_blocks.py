"""Matrix checks for the phase-polynomial intersection-block construction.

The script independently constructs (i) cyclic generator projectors and
(ii) branch-summed two-layer adaptive measurement effects.  It compares them
for the d=2,3,4 and d=6 examples discussed in the manuscript.
"""

from __future__ import annotations

import itertools
from collections.abc import Iterable

import numpy as np


Term = tuple[int, dict[int, int]]


def basis_index(values: Iterable[int], dimension: int) -> int:
    index = 0
    for value in values:
        index = index * dimension + value
    return index


def validate(
    dimension: int,
    parties: int,
    terms: list[Term],
    blocks: list[tuple[int, ...]],
    expected_gap: float,
) -> tuple[float, float]:
    omega = np.exp(2j * np.pi / dimension)
    basis = list(itertools.product(range(dimension), repeat=parties))
    hilbert_dimension = dimension**parties

    def polynomial(values: tuple[int, ...]) -> int:
        total = 0
        for coefficient, exponents in terms:
            monomial = coefficient
            for vertex, exponent in exponents.items():
                monomial *= values[vertex] ** exponent
            total += monomial
        return total % dimension

    diagonal = np.diag([omega ** polynomial(values) for values in basis])
    generator_projectors: list[np.ndarray] = []

    for vertex in range(parties):
        shift = np.zeros((hilbert_dimension, hilbert_dimension), dtype=complex)
        for values in basis:
            shifted = list(values)
            shifted[vertex] = (shifted[vertex] + 1) % dimension
            shift[
                basis_index(shifted, dimension),
                basis_index(values, dimension),
            ] = 1
        generator = diagonal @ shift @ diagonal.conj().T
        projector = sum(
            np.linalg.matrix_power(generator, power)
            for power in range(dimension)
        ) / dimension
        generator_projectors.append(projector)

    def algebraic_block(block: tuple[int, ...]) -> np.ndarray:
        effect = np.eye(hilbert_dimension, dtype=complex)
        for vertex in block:
            effect = effect @ generator_projectors[vertex]
        return effect

    def conditional_polynomial(
        vertex: int,
        value: int,
        complement_outcomes: dict[int, int],
    ) -> int:
        total = 0
        for coefficient, exponents in terms:
            if vertex not in exponents:
                continue
            monomial = coefficient * value ** exponents[vertex]
            for neighbour, exponent in exponents.items():
                if neighbour != vertex:
                    monomial *= complement_outcomes[neighbour] ** exponent
            total += monomial
        return total % dimension

    def measurement_block(block: tuple[int, ...]) -> np.ndarray:
        complement = tuple(v for v in range(parties) if v not in block)
        effect = np.zeros((hilbert_dimension, hilbert_dimension), dtype=complex)
        for outcomes in itertools.product(
            range(dimension), repeat=len(complement)
        ):
            measured = dict(zip(complement, outcomes, strict=True))
            passing_vector = np.zeros(hilbert_dimension, dtype=complex)
            for terminal_values in itertools.product(
                range(dimension), repeat=len(block)
            ):
                values = [0] * parties
                for vertex, outcome in measured.items():
                    values[vertex] = outcome
                phase = 0
                for vertex, value in zip(block, terminal_values, strict=True):
                    values[vertex] = value
                    phase += conditional_polynomial(vertex, value, measured)
                passing_vector[basis_index(values, dimension)] = (
                    omega**phase / dimension ** (len(block) / 2)
                )
            effect += np.outer(passing_vector, passing_vector.conj())
        return effect

    residuals = [
        np.linalg.norm(algebraic_block(block) - measurement_block(block), 2)
        for block in blocks
    ]
    verification_operator = sum(algebraic_block(block) for block in blocks) / len(
        blocks
    )
    eigenvalues = np.linalg.eigvalsh(verification_operator)[::-1]
    gap = float(np.real(1 - eigenvalues[1]))
    residual = float(max(residuals))

    assert residual < 1e-10
    assert abs(gap - expected_gap) < 1e-10
    return residual, gap


def main() -> None:
    four_party_terms: list[Term] = [
        (1, {0: 2, 2: 1}),
        (2, {1: 3, 2: 1, 3: 2}),
        (1, {2: 1, 3: 1}),
    ]
    for dimension in (2, 3, 4):
        residual, gap = validate(
            dimension,
            4,
            four_party_terms,
            [(0, 1), (2,), (3,)],
            1 / 3,
        )
        print(f"d={dimension}: residual={residual:.3e}, gap={gap:.12f}")

    residual, gap = validate(
        6,
        3,
        [(2, {0: 2, 2: 1}), (3, {1: 3, 2: 1})],
        [(0, 1), (2,)],
        1 / 2,
    )
    print(f"d=6: residual={residual:.3e}, gap={gap:.12f}")


if __name__ == "__main__":
    main()
