# Transfer and local-control checks

Run `python3 verify_resource_model.py` after the Li-Zhu branch files have been generated. It compares all 7368 exported branches at n=2,3,4 against an independent exact integer transfer calculation, verifies complete local-basis circuits and single-ray Givens rotations, and records example large-n single-record contractions.

See SI Note 11 for the recurrence, gate order, ideal continuous-control model, classical communication and precision assumptions. `compiled_three_bases.json` stores an initial diagonal operation followed by chronological two-level gates; every complex entry is `[real, imag]`. Fifteen rotations plus a diagonal phase layer are a generic complete-unitary upper bound, while a single accepted ray requires at most five rotations for either protocol. They are not hardware-optimal gate counts.
