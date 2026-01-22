SHELL := /bin/bash
.SHELLFLAGS := -ec
.ONESHELL:
.DEFAULT_GOAL := help

CONDA ?= mamba
CONDA_BUILD ?= conda
CONDA_CHANNELS ?= -c conda-forge -c bioconda -c defaults
ENV_FILE ?= environment.yml
ENV_NAME ?= scarab_cenv
TEST_ENV ?= scarab_test
PYTHON ?= python3
PY_VER ?= 3.11
PKG_NAME ?= scarab
RECIPE_LOCAL ?= conda-recipe-local
RECIPE_RELEASE ?= conda-recipe
IMAGE ?= quay.io/hallamlab/scarab
GIT_REF ?= $(shell git describe --tags --always --dirty 2>/dev/null || echo dev)
IMAGE_TAG ?= $(subst /,-,$(GIT_REF))
DOCKER ?= docker
DOCKER_SUDO ?= sudo
APPTAINER ?= apptainer
CONDA_BASE ?= $(shell $(CONDA) info --base 2>/dev/null)
CONDA_BLD_DIR ?= $(CONDA_BASE)/conda-bld
ANACONDA_USER ?=

.PHONY: help
help: ## Show available targets
	@awk 'BEGIN {FS = ":.*##"; print "Targets:"} /^[a-zA-Z0-9_-]+:.*##/ {printf "  %-24s %s\n", $$1, $$2}' $(MAKEFILE_LIST)

.PHONY: install-local
install-local: ## Install SCARAB into the active Python environment
	$(PYTHON) -m pip install .

.PHONY: install-editable
install-editable: ## Install SCARAB in editable mode into the active Python environment
	$(PYTHON) -m pip install -e .

.PHONY: env-create
env-create: ## Create a conda environment from environment.yml
	$(CONDA) env create -f $(ENV_FILE) -n $(ENV_NAME)

.PHONY: env-update
env-update: ## Update the conda environment from environment.yml
	$(CONDA) env update -f $(ENV_FILE) -n $(ENV_NAME) --prune

.PHONY: env-remove
env-remove: ## Remove the conda environment
	$(CONDA) env remove -n $(ENV_NAME)

.PHONY: conda-build-local
conda-build-local: ## Build from local sources using conda-recipe-local
	$(CONDA_BUILD) build --python $(PY_VER) $(RECIPE_LOCAL)

.PHONY: conda-build-release
conda-build-release: ## Build from the release recipe using conda-recipe
	$(CONDA_BUILD) build --python $(PY_VER) $(RECIPE_RELEASE)

.PHONY: conda-build-purge
conda-build-purge: ## Clear conda-build work and test directories
	$(CONDA_BUILD) build purge-all

.PHONY: conda-test-env
conda-test-env: ## Create a test env from the locally built package
	$(CONDA) create -n $(TEST_ENV) $(CONDA_CHANNELS) --use-local $(PKG_NAME)

.PHONY: conda-upload
conda-upload: ## Upload the latest noarch build to Anaconda Cloud (ANACONDA_USER required)
	if [ -z "$(ANACONDA_USER)" ]; then
		echo "ANACONDA_USER is not set"; exit 1;
	fi
	pkg=$$(ls -t $(CONDA_BLD_DIR)/noarch/$(PKG_NAME)-*.conda 2>/dev/null | head -1)
	if [ -z "$$pkg" ]; then
		pkg=$$(ls -t $(CONDA_BLD_DIR)/noarch/$(PKG_NAME)-*.tar.bz2 2>/dev/null | head -1)
	fi
	if [ -z "$$pkg" ]; then
		echo "No built package found in $(CONDA_BLD_DIR)/noarch"; exit 1;
	fi
	anaconda upload --user $(ANACONDA_USER) "$$pkg"

.PHONY: docker-build
docker-build: ## Build the Docker image from the current git ref
	$(DOCKER_SUDO) $(DOCKER) build --network=host \
		--build-arg git_branch=$(IMAGE_TAG) \
		-t $(IMAGE):$(IMAGE_TAG) .

.PHONY: docker-run
docker-run: ## Run the Docker image with the repo mounted at /cwd
	$(DOCKER_SUDO) $(DOCKER) run -it --network=host --rm \
		-v $(CURDIR):/cwd \
		$(IMAGE):$(IMAGE_TAG) bash

.PHONY: apptainer-build
apptainer-build: ## Build an Apptainer image from the Docker image
	$(APPTAINER) build $(PKG_NAME)-$(IMAGE_TAG).sif docker://$(IMAGE):$(IMAGE_TAG)

.PHONY: apptainer-run
apptainer-run: ## Run the Apptainer image
	$(APPTAINER) exec $(PKG_NAME)-$(IMAGE_TAG).sif bash
