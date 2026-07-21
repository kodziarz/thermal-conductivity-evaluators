Build with tests:

```shell
cmake -DBUILD_TESTS=ON -B ./build -S .
cmake --build ./build --config Release
```

Run the testsÖ
```
cd build
ctest
```