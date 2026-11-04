# Contributing

Thanks for your interest in improving the simulator.

## Development setup

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DRADAR_ENABLE_SANITIZERS=ON -DRADAR_WARNINGS_AS_ERRORS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Before opening a pull request

- Keep the code C++17, dependency-free and warning-free with `-Wall -Wextra -Wpedantic -Wshadow -Wconversion`.
- Add or update tests in `tests/` for every behaviour change. Tests compare against analytic
  expectations where possible (resolution, gains, Pfa), not against stored outputs.
- Run the static analysis used in CI:

  ```bash
  cppcheck --enable=warning,style,performance,portability --std=c++17 --inline-suppr \
    --suppressions-list=.cppcheck-suppressions --error-exitcode=1 -q -I include src app
  cmake -S . -B build-tidy -DCMAKE_CXX_COMPILER=clang++ -DRADAR_BUILD_TESTS=OFF
  run-clang-tidy -p build-tidy -quiet 'src/|app/'
  ```

- Keep results reproducible: new random draws must come from `mix_seed`-derived streams so the
  output does not depend on the OpenMP thread count.
- Update `README.md` and `CHANGELOG.md` (under *Unreleased*) for user-visible changes.

## Commit style

Small, focused commits with imperative subject lines, e.g. `Add Hann window option`.
