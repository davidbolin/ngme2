#!/usr/bin/env bash
set -euo pipefail

# sf calls a GDAL C++ entry point with std::unique_ptr in its signature.
# GDAL and sf must therefore use the same C++ standard library.
build_dir="$(mktemp -d /tmp/ngme2-gdal-clang.XXXXXX)"
curl --fail --location --retry 5 \
  https://github.com/OSGeo/gdal/releases/download/v3.12.4/gdal-3.12.4.tar.xz \
  --output "$build_dir/gdal.tar.xz"
printf '%s  %s\n' \
  813094498c17522ac42821a5ea1ea783d8326c0adf286cce86a949038bd09198 \
  "$build_dir/gdal.tar.xz" | sha256sum --check
tar -xf "$build_dir/gdal.tar.xz" -C "$build_dir"
cmake -S "$build_dir/gdal-3.12.4" -B "$build_dir/build" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=clang \
  -DCMAKE_CXX_COMPILER=clang++ \
  -DCMAKE_CXX_FLAGS=-stdlib=libc++ \
  -DCMAKE_INSTALL_PREFIX=/opt/gdal-clang \
  -DCMAKE_INSTALL_LIBDIR=lib \
  -DBUILD_APPS=OFF \
  -DBUILD_TESTING=OFF \
  -DGDAL_BUILD_OPTIONAL_DRIVERS=OFF \
  -DOGR_BUILD_OPTIONAL_DRIVERS=OFF
cmake --build "$build_dir/build" --parallel 2
cmake --install "$build_dir/build"
