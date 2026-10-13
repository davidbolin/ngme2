# Both compilers use the same Fedora 44 / R-devel base. The R-hub clang images
# use Ubuntu, so install Fedora's Clang and libc++ here instead.
ARG R_IMAGE=ghcr.io/r-hub/containers/gcc16@sha256:6238b814ab2be4feafbbb43e2e141ee8ecf647999753a9429ba90ef7a3a43264
FROM ${R_IMAGE}
ARG NGME_COMPILER=gcc

COPY Makevars.clang /opt/ngme2/Makevars.clang
COPY build-gdal-clang.sh /opt/ngme2/build-gdal-clang.sh
RUN if [ "$NGME_COMPILER" = clang ]; then \
      dnf install -y clang libcxx-devel libcxxabi-devel libomp-devel \
        cmake ninja-build proj-devel geos-devel sqlite-devel expat-devel \
        libtiff-devel libjpeg-turbo-devel libpng-devel zlib-devel && \
      bash /opt/ngme2/build-gdal-clang.sh && \
      cp /opt/ngme2/Makevars.clang /root/.R/Makevars; \
    elif [ "$NGME_COMPILER" != gcc ]; then \
      exit 2; \
    fi

ENV PATH=/opt/gdal-clang/bin:${PATH}
ENV LD_LIBRARY_PATH=/opt/gdal-clang/lib
ENV PKG_CONFIG_PATH=/opt/gdal-clang/lib/pkgconfig

COPY ngme2-source.tar.gz /tmp/ngme2-source.tar.gz
COPY install-deps.R /opt/ngme2/install-deps.R
# All R dependencies use the selected compiler and standard library. In the
# Clang image, sf also uses the matching GDAL installed above.
RUN Rscript /opt/ngme2/install-deps.R /tmp/ngme2-source.tar.gz
COPY check-toolchain.R /opt/ngme2/check-toolchain.R
RUN Rscript /opt/ngme2/check-toolchain.R
RUN mkdir -p /opt/ngme2-source && \
    tar -xzf /tmp/ngme2-source.tar.gz -C /opt/ngme2-source

COPY profile-tests.sh /usr/local/bin/profile-ngme2-tests
RUN chmod 0755 /usr/local/bin/profile-ngme2-tests

ENV NGME_TEST_SOURCE=/opt/ngme2-source/ngme2
ENV _R_CHECK_TESTS_NLINES_=0
ENV MAKEFLAGS=-j2
WORKDIR /check
CMD ["/usr/local/bin/profile-ngme2-tests"]
