# Structured fixed-measurement benchmarks

Run `python3 run_benchmark.py` for the physical matching checks, exact n=2 certificate and floating n=3 optimization. The script writes `benchmark_results.json` beside itself. Run `python3 exact_n2_certificate.py` to verify the rational certificate using the Python standard library alone.

The matching-only class, the CRT product-Pauli class, their analytic proofs and their precise scopes are defined in SI Note 10. In particular, the 1/6 matching bound is not a bound on all nonadaptive local measurements. The n=3 search in `pauli_search.py` uses NumPy and exports feasible numerical primal and dual candidates; it is not an interval-certified SDP solver.
