# Regenerate package data for documentation examples.
# This is a maintainer-only script; qcc is not a runtime package dependency.

data("pistonrings", package = "qcc")

if (!dir.exists("data")) {
  dir.create("data")
}

save(pistonrings, file = "data/pistonrings.rda", compress = "xz")
