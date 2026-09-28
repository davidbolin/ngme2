FROM ghcr.io/r-hub/containers/gcc16@sha256:6aa66e5f02f9cbca5eba368c532b612dcdbbdd817f16cca48f9227f7d1a7cfb6

COPY ngme2_1.0.0.tar.gz /tmp/ngme2_1.0.0.tar.gz

# Match a source-package check: install declared dependencies, then compile
# the exact tarball that was submitted to CRAN.
RUN Rscript -e 'pak::pkg_install("deps::/tmp/ngme2_1.0.0.tar.gz", dependencies = TRUE)'
RUN R CMD INSTALL /tmp/ngme2_1.0.0.tar.gz
RUN mkdir -p /opt/ngme2-source && tar -xzf /tmp/ngme2_1.0.0.tar.gz -C /opt/ngme2-source

COPY profile-tests.sh /usr/local/bin/profile-ngme2-tests
RUN chmod 0755 /usr/local/bin/profile-ngme2-tests

ENV NGME_TEST_SOURCE=/opt/ngme2-source/ngme2
WORKDIR /check
CMD ["/usr/local/bin/profile-ngme2-tests"]
