# SeaExplorer Glider Data Processing - Scientific Documentation

This documentation provides a comprehensive guide to processing and quality controlling SeaExplorer glider data using PyGlider, following best practices for oceanographic data management.

## Documentation Structure

### 1. Main Article
**File:** `best_practices_article.md`

Formatted as a scientific article suitable for submission to *Frontiers in Marine Science* (Best Practices section), including:
- Abstract
- Introduction
- Methods (detailed processing pipeline)
- Results (validation and quality metrics)
- Discussion
- Data Availability Statement
- Code Availability
- References

### 2. Technical Manual
**File:** `technical_manual.md`

Detailed technical documentation covering:
- System requirements
- Installation procedures
- Complete pipeline workflow
- Configuration guidelines
- Troubleshooting

### 3. Quality Control Guide
**File:** `quality_control_guide.md`

Comprehensive QC documentation including:
- QC flag definitions
- Range specifications
- Sensor specifications
- Geographic constraints (LAND_QC)
- Validation procedures

### 4. Unit Conversions Reference
**File:** `unit_conversions.md`

Detailed documentation of all unit conversions applied:
- Turbidity (backscatter → NTU)
- Conductivity (mS/cm → S/m)
- Dissolved oxygen (mmol/L → µmol/kg)
- Salinity recalculation (TEOS-10)
- Density recalculation (TEOS-10)

### 5. Code Examples
**Directory:** `examples/`

Practical code examples and use cases:
- Basic pipeline execution
- Custom QC implementation
- Data visualization
- Integration with other tools

## Target Audience

- Marine scientists working with autonomous underwater vehicles
- Data managers handling glider datasets
- Oceanographers requiring standardized data processing
- Researchers preparing data for publication in data repositories

## Citation

If you use this processing pipeline in your research, please cite:

```
Torelli, B. (2025). Best Practices for SeaExplorer Glider Data Processing 
and Quality Control Using PyGlider. GitHub repository. 
https://github.com/BennyTorelli/Seaexplorer_pyglider
```

## Contributing

This documentation is maintained as part of the SeaExplorer_PyGlider project. 
Contributions, corrections, and suggestions are welcome via GitHub issues or pull requests.

## License

This documentation is released under Creative Commons Attribution 4.0 International (CC BY 4.0).

## Contact

For questions or collaboration inquiries regarding this documentation, please open an issue on the GitHub repository.

---

**Repository:** https://github.com/BennyTorelli/Seaexplorer_pyglider  
**Last Updated:** October 2025
