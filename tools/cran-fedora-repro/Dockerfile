# Both compilers use the same Fedora 44 / R-devel base. The R-hub clang images
# use Ubuntu, so install Fedora's Clang and libc++ here instead.
ARG R_IMAGE=ghcr.io/r-hub/containers/gcc16@sha256:6238b814ab2be4feafbbb43e2e141ee8ecf647999753a9429ba90ef7a3a43264
FROM ${R_IMAGE}
ARG NGME_COMPILER=gcc

COPY Makevars.clang /opt/ngme2/Makevars.clang
RUN if [ "$NGME_COMPILER" = clang ]; then \
      dnf install -y clang libcxx-devel libcxxabi-devel libomp-devel && \
      cp /opt/ngme2/Makevars.clang /root/.R/Makevars; \
    elif [ "$NGME_COMPILER" != gcc ]; then \
      exit 2; \
    fi

COPY ngme2-source.tar.gz /tmp/ngme2-source.tar.gz
COPY install-deps.R /opt/ngme2/install-deps.R
RUN Rscript /opt/ngme2/install-deps.R /tmp/ngme2-source.tar.gz
RUN mkdir -p /opt/ngme2-source && \
    tar -xzf /tmp/ngme2-source.tar.gz -C /opt/ngme2-source

COPY profile-tests.sh /usr/local/bin/profile-ngme2-tests
RUN chmod 0755 /usr/local/bin/profile-ngme2-tests

ENV NGME_TEST_SOURCE=/opt/ngme2-source/ngme2
ENV _R_CHECK_TESTS_NLINES_=0
ENV MAKEFLAGS=-j2
WORKDIR /check
CMD ["/usr/local/bin/profile-ngme2-tests"]
