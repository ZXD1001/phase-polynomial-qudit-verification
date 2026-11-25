# Optimized comparison certificates

Run through `../run_revision_checks.py` to preserve these sources and inputs. It copies this directory to a fresh run and executes the five preserved reference checks plus `verify_n4_primal.py`.

- `exact_certificate.py`: exact n=2 CRT reference optimum 3/11 in the prescribed 144-setting class.
- `certify_lz_n4_upper.py --verify-only`: exact rational dual upper `197537/500000` for all 32 terminal/Z-F-front tests.
- `certify_fullmenu_half.py --verify-only`: exact half gap using 23 settings whose first-layer PVMs are in the four-PVM library but terminal rays are unrestricted. This is not a fixed-complete-terminal-PVM width-one protocol.
- `check_identities.py`: native n=3 operator residuals and uniform regression reconstruction; finite exact binary-rank/compatibility checks support the separately stated analytic arguments.
- `count_certificate_rays.py`: ray/PVM counts for the supplied implementation only.
- `verify_n4_primal.py`: floating re-evaluation of the supplied n=4 optimized weights, without an optimizer. Numerical tolerance 1e-12; not an interval certificate.

`lz_optimize.py` is preserved as the reconstruction dependency used by the identity/primal checks. Its optimization entry point is not run. No uniform table is interpreted as an optimized method limit.
