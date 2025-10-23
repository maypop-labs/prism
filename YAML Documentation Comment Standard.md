# YAML Documentation Comment Standard v1.0

## Overview

This document defines a comprehensive tag-based comment system for YAML configuration files that enables automatic documentation generation and provides detailed context for configuration parameters. The system uses structured comment tags that can be parsed by documentation generators while maintaining human readability.

## Core Principles

1. **Machine Parseable**: All documentation uses `@tag:` format for easy extraction
2. **Human Readable**: Comments remain clear and useful when read directly
3. **Practical Guidance**: Provides actionable recommendations with rationales
4. **Type-Safe**: Documents expected data types and constraints
5. **Maintainable**: Facilitates long-term project maintenance and onboarding
6. **Honest Documentation**: Only claims enforcement when actually implemented
7. **Semantic Clarity**: Distinct tags for different types of information

## Tag Specification

### Section-Level Tags

These tags define and describe major configuration sections:

#### `@doc: [section_name]`
- **Purpose**: Marks the beginning of a new configuration section
- **Format**: Single word or underscore_separated identifier
- **Usage**: Place before the first parameter in a logical grouping
- **Example**: `@doc: File_Paths`

#### `@info: [description]`
- **Purpose**: Describes what the entire section configures
- **Format**: Free-form text describing section purpose
- **Usage**: Follows `@doc:` tag immediately
- **Example**: `@info: Core directory paths for data input/output operations`

### Parameter-Level Tags

These tags document individual configuration parameters:

#### `@param: [parameter_name]`
- **Purpose**: Identifies the specific parameter being documented
- **Format**: Exact parameter name as it appears in YAML
- **Usage**: Must match the YAML key exactly
- **Example**: `@param: seuratMaxPercentMtDNA`

#### `@type: [data_type]`
- **Purpose**: Specifies the expected data type and format
- **Format**: Use standard types with clarifications in parentheses
- **Valid Types**:
  - `string` - Text values
  - `string (path)` - File/directory paths
  - `integer` - Whole numbers
  - `float` - Decimal numbers
  - `boolean` - true/false values
  - `list` - Array of values
  - `list[string]` - Array of specific type
- **Example**: `@type: integer`

#### `@description: [comprehensive_description]`
- **Purpose**: Provides comprehensive explanation of parameter function, including warnings and constraints
- **Format**: Complete sentences with optional multi-line continuation
- **Usage**: Should include all relevant information (former @note content integrated here)
- **Multi-line Format**: Continuation lines start with `# ` (no extra indentation)
- **Example**: 
  ```yaml
  # @description: Maximum percentage of mitochondrial genes allowed per cell during quality control.
  # Cells exceeding this threshold are filtered out as potentially dying/damaged.
  ```

### Guidance Tags

These tags provide practical recommendations and insights:

#### `@suggested: [recommended_values]`
- **Purpose**: Provides recommended values with options for different use cases
- **Format**: List format with conservative, standard, and/or aggressive options
- **Usage**: Gives users practical starting points based on common scenarios
- **Example**: `@suggested: [2, 4, 8]` or `@suggested: [0.1, 0.05, 0.01]`

#### `@rationale: [explanation_of_suggestions]`
- **Purpose**: Explains the reasoning behind suggested values and trade-offs
- **Format**: Educational explanation of when to use different values
- **Usage**: Helps users make informed decisions about parameter choices
- **Example**: `@rationale: 2 cores for basic processing, 4 cores for balanced performance, 8 cores for compute-intensive analyses. Values above system core count provide no benefit.`

### Optional Parameter Tags

These tags provide additional context when needed:

#### `@default: [default_value]`
- **Purpose**: Documents default values when parameter is optional
- **Format**: The actual default value used if parameter is omitted
- **Usage**: Clarifies behavior when parameter is not specified
- **Example**: `@default: 1000`

#### `@range: [constraints]`
- **Purpose**: Documents enforced ranges or constraints in the code
- **Format**: Range notation or descriptive text
- **Usage**: **Only use when ranges are actually validated in the code**
- **Example**: `@range: 0.0-1.0` (only if code actually checks this)
- **Important**: Do not use for suggested ranges - use `@suggested` and `@rationale` instead

## Comment Structure Template

```yaml
# @doc: [Section_Name]
# @info: [What this section configures]

# @param: parameter_name
# @type: [data_type]
# @description: [Detailed description of what this parameter does]
# [Optional continuation lines for comprehensive information]
# @suggested: [recommended_values]
# @rationale: [Why these values are recommended and when to use each]
# @default: [default_if_optional]
# @range: [only_if_actually_enforced]
parameter_name: actual_value
```

## Multi-Line Descriptions

The `@description:` tag supports multi-line text for comprehensive documentation:

### Format
```yaml
# @description: Primary description of the parameter function.
# Additional details, warnings, or constraints.
# Special considerations or dependencies.
```

### Guidelines
- **No Extra Indentation**: Continuation lines use only `# ` (hash + single space)
- **Natural Flow**: Information should flow logically from general to specific
- **Comprehensive**: Include all relevant information from former separate notes
- **Maintainable**: Easy to edit without worrying about alignment

### Parsing Logic
1. When encountering `@description:`, capture that line's content
2. Continue reading subsequent lines that start with `# ` and don't contain a new `@tag:`
3. Concatenate all lines as the complete description

## Usage Guidelines

### Tag Ordering
1. Section tags (`@doc:`, `@info:`) appear first for each section
2. Parameter tags appear immediately before each parameter
3. Parameter tags follow this order:
   - `@param:` (required)
   - `@type:` (required)
   - `@description:` (required, may be multi-line)
   - `@suggested:` (recommended)
   - `@rationale:` (when @suggested is used)
   - `@default:` (when applicable)
   - `@range:` (only when enforced)

### Content Guidelines

#### For @description
- **Be Comprehensive**: Include all relevant information in one place
- **Include Constraints**: Document important limitations or requirements
- **Warn About Dependencies**: Note when parameters depend on external files or other parameters
- **Use Natural Language**: Write as if explaining to a colleague

#### For @suggested and @rationale
- **Provide Options**: Give 2-3 values covering conservative, standard, and aggressive use cases
- **Explain Trade-offs**: Help users understand the implications of each choice
- **Be Educational**: Teach users about the underlying concepts
- **Consider Context**: Account for different computational resources and use cases
- **Be Honest**: Don't claim enforcement that doesn't exist

#### General Guidelines
- **Be Specific**: Use precise language that explains both what and why
- **Include Context**: Explain how parameters affect pipeline behavior
- **Update Regularly**: Keep documentation current when code changes
- **Test Examples**: Verify that suggested values actually work

### Section Organization
Group related parameters logically:
- **File_Paths**: File and directory locations
- **Experiment_Design**: Sample metadata and experimental setup
- **System_Resources**: System resources and runtime settings
- **Quality_Control**: Cell and gene filtering parameters
- **Analysis_Stage_Name**: Parameters specific to each analysis stage
- **Visualization**: Plotting and output formatting

## Example Implementation

```yaml
# @doc: System_Resources
# @info: System resource allocation and computation environment settings

# @param: cores
# @type: integer
# @description: Number of CPU cores to use for parallel processing operations.
# @suggested: [2, 4, 8]
# @rationale: 2 cores for basic processing, 4 cores for balanced performance, 8 cores for compute-intensive analyses. Values above system core count provide no benefit.
# @default: 1
cores: 4

# @param: seuratMaxPercentMtDNA
# @type: float
# @description: Maximum percentage of mitochondrial genes allowed per cell during quality control.
# Cells exceeding this threshold are filtered out as potentially dying/damaged.
# Must be set before running quality control steps.
# @suggested: [15.0, 20.0, 25.0]
# @rationale: 15% for stringent filtering (removes more potentially damaged cells), 20% for standard filtering (balances quality and cell retention), 25% for permissive filtering (retains more cells but may include damaged ones).
# @default: 20.0
seuratMaxPercentMtDNA: 20.0

# @param: rootPath
# @type: string (path)
# @description: Base directory containing cellranger_counts folder and where all results are written.
# Must contain cellranger_counts subfolder with 10X Genomics output.
rootPath: "E:/datasets/omics/skin/"
```

## Documentation Generation

This standard enables automatic generation of:
1. **Parameter Reference**: Complete list with descriptions and recommendations
2. **Decision Guide**: Suggested values with rationales for informed choices
3. **Type Reference**: Data types and validation rules
4. **Configuration Templates**: Starting configurations for different use cases
5. **Educational Resources**: Understanding the science behind parameter choices

## Best Practices

### For @description
- **Start General**: Begin with the primary function of the parameter
- **Add Specifics**: Include constraints, warnings, and dependencies
- **Use Plain Language**: Write for both experts and newcomers
- **Keep Current**: Update when parameter behavior changes

### For @suggested and @rationale
- **Test Your Suggestions**: Only recommend values you've actually used successfully
- **Explain Edge Cases**: Note when extreme values might be appropriate
- **Consider Resources**: Account for different computational environments
- **Update Based on Experience**: Refine recommendations as you learn more

### General Best Practices
- **Start with Core Tags**: Begin with @param, @type, @description, then add guidance
- **Stay Current**: Update documentation when code behavior changes
- **Be Complete**: Document every parameter, even "obvious" ones
- **Test Documentation**: Verify examples and suggestions actually work
- **Review Periodically**: Ensure documentation matches current understanding

### Common Patterns

#### For Resource Parameters
```yaml
# @suggested: [conservative, standard, intensive]
# @rationale: Conservative for limited resources, standard for typical usage, intensive for high-performance systems.
```

#### For Threshold Parameters
```yaml
# @suggested: [strict, moderate, permissive]
# @rationale: Strict for high-quality results, moderate for balanced approach, permissive for exploratory analysis.
```

#### For File Paths
```yaml
# @description: [Path description].
# Must exist and contain required subfolders.
# @required: true
```

This standard transforms configuration files from simple parameter lists into educational resources that help users make informed decisions about their analyses while maintaining clean, unambiguous documentation structure.