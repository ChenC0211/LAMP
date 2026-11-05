# Formal FirstJet Hydro App fixture

These seven dictionaries preserve the original uniform 32 x 32 static-drop input used by `hydro_app_contract.py`: Hermite-GALS, staggered Hydro, cell-average momentum density, balanced CSF, and exact zero viscosity.

They were copied byte for byte from the pre-VOF-migration `cases/hydro_static_drop/system` backup on 2026-10-17. `provenance.json` records the exact backup location and SHA256 of every original dictionary. Keep numerical values and acceptance checks independent of later production case changes.

The FirstJet-specific `fresh()` and index-lock contracts copy `system/`. The four standard-case preflights continue to copy the actual `cases/<name>/system` inputs. This fixture is test input, not a production case or a VOF validation window.
