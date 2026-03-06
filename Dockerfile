# syntax=docker/dockerfile:1

FROM vastai/pytorch:cuda-12.4.1-auto

ENV DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    CUDA_HOME=/usr/local/cuda \
    PIP_EXTRA_INDEX_URL=https://download.pytorch.org/whl/cu121 \
    HF_HOME=/cache
    
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    software-properties-common \
    curl \
    && add-apt-repository -y ppa:deadsnakes/ppa \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
    python3.11 \
    python3.11-dev \
    python3.11-venv \
    build-essential \
    git \
    cmake \
    pkg-config \
    libffi-dev \
    libssl-dev \
    libxml2-dev \
    libxslt-dev \
    libgl1 \
    libhdf5-dev \
    libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.11 1 && \
    update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1 && \
    python3.11 -m pip install --upgrade --ignore-installed setuptools setuptools_scm wheel

WORKDIR /app

COPY . /app

RUN pip install --no-cache-dir -e /app

ARG DOWNLOAD_WEIGHTS=false
RUN mkdir -p "${HF_HOME}" && \
    if [ "${DOWNLOAD_WEIGHTS}" = "true" ]; then \
        boltzgen download all --cache "${HF_HOME}" --force_download; \
    fi

ARG USERNAME=boltzgen
ARG USER_UID=1000
ARG USER_GID=1000

RUN if ! getent group "${USER_GID}" >/dev/null; then groupadd --gid "${USER_GID}" "${USERNAME}"; fi && \
    if ! id -u "${USERNAME}" >/dev/null 2>&1; then \
        if getent passwd "${USER_UID}" >/dev/null; then \
            useradd --gid "${USER_GID}" --create-home --shell /bin/bash "${USERNAME}"; \
        else \
            useradd --uid "${USER_UID}" --gid "${USER_GID}" --create-home --shell /bin/bash "${USERNAME}"; \
        fi; \
    fi

RUN mkdir -p "${HF_HOME}" && chown -R ${USER_UID}:${USER_GID} "${HF_HOME}"

USER ${USERNAME}
WORKDIR /workspace
