SHELL := /bin/bash

.PHONY: help test-quick test-unit test-all lint smoke smoke-mock smoke-real smoke-poor smoke-all dev-setup

help: ## List available targets
	@echo "Available targets:"
	@awk 'BEGIN {FS = ":.*## "} /^[a-zA-Z0-9_.-]+:.*## / {printf "  %-16s %s\n", $$1, $$2}' $(MAKEFILE_LIST)

# =============================================================================
# Unit & Integration Tests
# =============================================================================

test-quick: ## Run quick tests (unit tests only, no data generation)
	bash tests/run_tests.sh --quick

test-unit: ## Run unit tests only
	bash tests/run_tests.sh --unit-only

test-all: ## Run full test suite
	bash tests/run_tests.sh

# =============================================================================
# Smoke Tests (End-to-End)
# =============================================================================

smoke: smoke-mock ## Run smoke test with mock tools (default, fast)

smoke-mock: ## Run smoke test with mock tools (~10s, no deps required)
	@echo "Running mock smoke test..."
	bash tests/smoke/run_smoke.sh --mock --profile good

smoke-real: ## Run smoke test with real tools (requires bwa, samtools, etc.)
	@echo "Running real smoke test..."
	bash tests/smoke/run_smoke.sh --real --profile good --verbose

smoke-poor: ## Run smoke test with poor quality data (edge case testing)
	@echo "Running poor-quality smoke test..."
	bash tests/smoke/run_smoke.sh --mock --profile poor

smoke-all: ## Run smoke tests for all quality profiles
	@echo "Running all smoke tests..."
	@for profile in good poor adapter mixed; do \
		echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"; \
		echo "Testing profile: $$profile"; \
		echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"; \
		bash tests/smoke/run_smoke.sh --mock --profile $$profile || exit 1; \
	done
	@echo "All smoke profiles passed!"

smoke-keep: ## Run smoke test and keep output artifacts for inspection
	bash tests/smoke/run_smoke.sh --mock --profile good --keep --verbose

# =============================================================================
# Development
# =============================================================================

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

dev-setup: ## Install development dependencies (pre-commit)
	pip3 install pre-commit
	pre-commit install

# =============================================================================
# CI Targets
# =============================================================================

ci-quick: test-unit smoke-mock ## Quick CI validation (unit tests + mock smoke)
	@echo "CI quick validation passed!"

ci-full: test-all smoke-all ## Full CI validation (all tests + all smoke profiles)
	@echo "CI full validation passed!"
