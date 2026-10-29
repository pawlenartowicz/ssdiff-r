test_that("need_pkg names the function, the package and install_backends()", {
  pkg <- "ssdiffNoSuchPackage"
  err <- tryCatch(need_pkg(pkg, "fit_pls()"), error = function(e) e)
  expect_s3_class(err, "error")
  msg <- conditionMessage(err)
  expect_match(msg, "fit_pls() needs the 'ssdiffNoSuchPackage' package, which", fixed = TRUE)
  expect_match(msg, "ssdiff::install_backends()", fixed = TRUE)
  expect_match(msg, "install.packages(\"ssdiffNoSuchPackage\", repos = c(\"https://pawlenartowicz.r-universe.dev\"",
               fixed = TRUE)
  expect_null(conditionCall(err))
})

test_that("need_pkg puts the reason in the message", {
  expect_error(need_pkg("ssdiffNoSuchPackage", "fit_regression()", "a formula with random effects"),
               "fit_regression() needs the 'ssdiffNoSuchPackage' package (for a formula with random effects), which",
               fixed = TRUE)
})

test_that("need_pkg is silent when the package loads", {
  expect_silent(need_pkg("stats", "fit_regression()"))
})

test_that("need_pkg stops on a version older than DESCRIPTION's minimum and names both", {
  local_mocked_bindings(.min_version = function(pkg) "999.0")
  err <- tryCatch(need_pkg("stats", "fit_regression()"), error = function(e) conditionMessage(e))
  expect_match(err, sprintf("fit_regression() needs 'stats' 999.0 or later; version %s is installed.",
                            utils::packageVersion("stats")), fixed = TRUE)
  expect_match(err, "ssdiff::install_backends()", fixed = TRUE)
})

test_that("the minimum versions come from DESCRIPTION's Suggests", {
  suggests <- utils::packageDescription("ssdiff", fields = "Suggests")
  for (pkg in c("plskit", "fastglmm")) {
    expect_match(.min_version(pkg), "^[0-9]+(\\.[0-9]+)+$")
    expect_match(suggests, sprintf("%s (>= %s)", pkg, .min_version(pkg)), fixed = TRUE)
  }
  expect_null(.min_version("stats"))
})
