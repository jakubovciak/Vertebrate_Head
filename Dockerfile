FROM bioconductor/bioconductor_docker:RELEASE_3_15

LABEL name="Markos 2024 Transitions" \
    version=0.1 \
    maintainer="kubovcij@img.cas.cz" \
    description="Data, code and software dependencies for Markos 2023 transitions analysis" \
    license="MIT"


ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

RUN apt-get update && apt-get install -y \
    texlive-latex-extra \
    texlive-fonts-recommended \
    dvipng \
    pandoc

RUN R -e "BiocManager::install('rstudio/renv@0.15.5')"

COPY Transitions/software/renv.lock /

RUN Rscript -e "print(BiocManager::repositories())"

RUN Rscript -e "\
    options(Ncpus = 4, pkgType = 'binary');\
    renv::consent(TRUE);\
    renv::restore(lockfile = 'renv.lock', prompt = FALSE, repos = BiocManager::repositories());\
    "
	
COPY Transitions/software/environment.yml /

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_23.11.0-2-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-py39_23.11.0-2-Linux-x86_64.sh -b \
    && rm -f Miniconda3-py39_23.11.0-2-Linux-x86_64.sh
	
RUN conda config --set ssl_verify false

RUN conda env create -y -f environment.yml

RUN conda init

CMD ["/init"]
