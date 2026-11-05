## Summary

<!-- What does this change and why? -->

## Verification

<!-- Tests added or updated, analytic expectation checked, scenarios run. -->

## Checklist

- [ ] `ctest` passes in a sanitizer build (`-DRADAR_ENABLE_SANITIZERS=ON`)
- [ ] No new warnings with `-DRADAR_WARNINGS_AS_ERRORS=ON` (GCC and Clang)
- [ ] cppcheck and clang-tidy are clean
- [ ] README and CHANGELOG updated for user-visible changes
