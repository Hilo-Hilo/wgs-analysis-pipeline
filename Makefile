SHELL := /bin/bash

.PHONY: help test-quick test-unit test-all lint

help: ## List available targets
	@echo "Available targets:"
	@awk 'BEGIN {FS = ":.*## "} /^[a-zA-Z0-9_.-]+:.*## / {printf "  %-12s %s\n", $$1, $$2}' $(MAKEFILE_LIST)

test-quick: ## Run quick tests
	bash tests/run_tests.sh --quick

test-unit: ## Run unit tests only
	bash tests/run_tests.sh --unit-only

test-all: ## Run full test suite
	bash tests/run_tests.sh

lint: ## Run shellcheck on all .sh files (if installed)
	@if command -v shellcheck >/dev/null 2>&1; then \
		echo "Running shellcheck..."; \
		scripts="$$(find . -type f -name '*.sh' -not -path './.git/*')"; \
		if [[ -z "$$scripts" ]]; then \
			echo "No shell scripts found."; \
		else \
			shellcheck -S error -e SC2034,SC2155,SC2086 $$scripts; \
		fi; \
	else \
		echo "shellcheck is not installed; skipping lint."; \
	fi
