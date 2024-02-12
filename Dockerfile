FROM bioconductor/bioconductor_docker:RELEASE_3_15

LABEL name="Markos 2024 Transitions" \
    version=0.1 \
    maintainer="kubovcij@img.cas.cz" \
    description="Data, code and software dependencies for Markos 2023 transitions analysis" \
    license="MIT"

RUN R -e "BiocManager::install('rstudio/renv@0.15.5')"

COPY Transitions/software/renv.lock /

RUN Rscript -e "print(BiocManager::repositories())"

RUN Rscript -e "\
    options(Ncpus = 4, pkgType = 'binary');\
    renv::consent(TRUE);\
    renv::restore(lockfile = 'renv.lock', prompt = FALSE, repos = BiocManager::repositories());\
    "

CMD ["/init"]
