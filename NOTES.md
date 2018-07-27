# This file contains a few notes on what work has to be done in this repo

## Urgent work
- Update the doc to be publishable
- Choose a license

## New features
- Test current implementation for supported types combinations
- Make sure everything that is advertised in the manual is implemented
- Implement a test suite to look for bugs

## Design changes

## Minor changes
- Change the types used. This means:
  - Remove usage of int in favor of `std::int64_t` or `std::uint64_t`.
  - Change `std::int64_t` to `std::uint64_t` in some places.
