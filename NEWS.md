# pboost4sar 0.6.2 (2026-04-22)

* add `forward-sar.R`.


# pboost4sar 0.6.1 (2026-01-15)

* add attrubite `minusloglik` in `get_rho()`, and improve `psar()` acccordingly.



# pboost4sar 0.6.0 (2026-01-14)

* rename the package as `pboost4sar`.

# pvs4sar 0.6.0 (2026-01-13)

* improve the computational efficiency of adaptive lasso for SAR model, including:
    - use Woodbury matrix identity for ridge estimation
    - add arguments `beta.ini` and `rho.ini` for grid search
    - remove the shared constant `log(det(A.rho))` in BIC



# pvs4sar 0.5.2 (2025-12-25)

* add `flag` into results of all methods.


# pvs4sar 0.5.1 (2025-12-25)

* update `simu-sar-data.R` according to the manuscript.




# pvs4sar 0.5.0 (2025-12-23)

* mv `R/plagsarlm.R` to `R/pvs4sar.R`.

* mv `R/sam-adaptivelasso.R` to `R/penalized-sar-adaptivelasso.R`.

* add `R/penalized-sar.R`, `R/penalized-sar-lasso.R`, `R/penalized-sar-scad.R`.


# pvs4sar 0.4.0 (2025-12-16)

* rename the package name as "pvs4sar" (from "pboostsam").

* eval argument `data` at first in `plagsarlm()`.

* rename `simu_sam_data_*()` as `simu_sar_data_*()` and the releated filenames and examples.



# pboostsam 0.3.2 (2025-12-15)

* `R/sam-adaptivelasso.R`: add argument `criterion = "BIC" / "EBIC"` in `tune_sam_adaptivelasso()`.

* improve stability of some `stopifnot`.




# pboostsam 0.3.1 (2025-12-14)

* add argument `rho` to `tune_sam_adaptivelasso()`.

* polish documentation.



# pboostsam 0.3.0 (2025-12-08)

* rename `pboostsam.R` as `plagsarlm.R`.

* add `plagsarlm2()` in `plagsarlm.R`.

* add `sam-adaptivelasso.R`.

* add `get-rho.R`.



# pboostsam 0.2.2 (2025-11-25)

* polish the generaction of rook weight matrix in `simu_sam_data_rook()` in `simu-data-data.R`.


# pboostsam 0.2.1 (2025-11-21)

* add `simu_sam_data_case()` in `simu-sam-data.R`.



# pboostsam 0.2.0 (2025-11-21)

* fix the generation of rook matrix in `simu-sam-data.R`.



# pboostsam 0.1.0 (2025-11-20)

* Initialization.
