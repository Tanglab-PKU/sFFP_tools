# QC Model for Cross-Linking Mass Spectrometry Analysis

This repository contains machine learning models and scripts for analyzing cross-linking mass spectrometry (XL-MS) data, specifically designed to work with PRIDE datasets and pLink3 results.

## Overview

### AI_model
Contains a Random Forest classifier trained on PRIDE data (PXD048452) for photo-cross-link identification. The model incorporates multiple features including:
- Intensity distributions (mean and standard error) of backbone (b/y) and side-chain (sc/sz) fragments
- Counts of characteristic fragment ions
- Peptide-spectrum match (PSM) scores

**Key Model Specifications:**
- Random Forest classifier with 500 trees
- Maximum depth = 4
- Feature importance analysis revealed PSM counts, backbone b/y ion coverage, and sz fragment ion intensity as the most discriminative features

### simplified_model_workflow
A simplified version that focuses on the core features of the Random Forest model. This version is recommended for most use cases and has been tested with:
- pLink3.0.16
- MeroX2.1.0
- xiSEARCH1.8.7

## Installation & Requirements

- **Python Version:** 3.7 or higher

## Usage

### Basic Workflow

1. Download the corresponding PRIDE dataset
2. Run the scripts in numerical order
