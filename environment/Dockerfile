# hash:sha256:ceaac4e158e1ec4000763fa4a7125cf9f33991b0951ec4acfd2468ee28546255
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/c1-vscode:4.20.0-pylint2023.10.1-copilot1.143.601

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-ask-pass /

RUN mamba install -y -c conda-forge \
    python==3.11 \
    ipykernel==6.29.5 \
    && conda clean -ya

# system libs for zarr / scikit-image / scipy (OpenMP)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgl1 libglib2.0-0 libgomp1
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

# postInstall clones + installs 2p2xenium (the xenium_autocoreg package). Its pyproject
# pins the runtime deps (numpy/pandas/scipy/scikit-image/tifffile/zarr/matplotlib/tqdm) --
# no extra pip installs needed here.
COPY postInstall /
RUN --mount=type=secret,id=secrets . /run/secrets/secrets \
    && /postInstall
