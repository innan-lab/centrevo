# Contributing to Centrevo

Thank you for your interest in contributing to Centrevo! This document provides guidelines and instructions for contributing.

## 🌟 How to Contribute

There are many ways to contribute to Centrevo:

- 🐛 Report bugs
- 💡 Suggest new features or enhancements
- 📝 Improve documentation
- 🧪 Add tests
- 🔧 Fix bugs or implement features
- 📊 Add benchmarks
- 🎨 Improve code quality

## 🚀 Getting Started

### Prerequisites

- Rust 1.70 or later
- Git
- (Optional) Python 3.8+ and maturin for Python bindings

### Setup Development Environment

1. **Fork and clone the repository:**
   ```bash
   git clone https://github.com/YOUR_USERNAME/centrevo.git
   cd centrevo
   ```

2. **Create a branch:**
   ```bash
   git checkout -b feature/your-feature-name
   ```

3. **Build the project:**
   ```bash
   cargo build
   ```

4. **Run tests:**
   ```bash
   cargo test
   ```

## 📋 Pull Request Process

1. **Update your branch:**
   ```bash
   git fetch origin
   git rebase origin/main
   ```

2. **Make your changes:**
   - Write clear, concise commit messages
   - Follow the code style guidelines (see below)
   - Add tests for new functionality
   - Update documentation as needed

3. **Test your changes:**
   ```bash
   cargo test
   cargo clippy
   cargo fmt --check
   ```

4. **Push your branch:**
   ```bash
   git push origin feature/your-feature-name
   ```

5. **Create a Pull Request:**
   - Go to the GitHub repository
   - Click "New Pull Request"
   - Select your branch
   - Fill in the PR template
   - Request review

6. **Address feedback:**
   - Respond to code review comments
   - Make requested changes
   - Push updates to your branch

## 📝 Code Style Guidelines

### Rust Code

- **Follow Rust conventions:**
  - Use `snake_case` for functions and variables
  - Use `PascalCase` for types and traits
  - Use `SCREAMING_SNAKE_CASE` for constants

- **Code formatting:**
  ```bash
  cargo fmt
  ```

- **Linting:**
  ```bash
  cargo clippy
  ```
  All code should pass clippy with no warnings.

- **Documentation:**
  - Add doc comments (`///`) for public APIs
  - Include examples in doc comments when helpful
  - Use `//!` for module-level documentation

- **Error handling:**
  - Use `Result` for fallible operations
  - Create custom error types when appropriate
  - Provide helpful error messages

### Testing

- **Write tests for new features:**
  ```rust
  #[cfg(test)]
  mod tests {
      use super::*;

      #[test]
      fn test_your_feature() {
          // Test code
      }
  }
  ```

- **Test naming:**
  - Use descriptive names: `test_sequence_mutation_changes_bases`
  - Follow pattern: `test_<component>_<scenario>_<expected_result>`

- **Test coverage:**
  - Aim for high coverage of new code
  - Test edge cases and error conditions
  - Include integration tests when appropriate

### Benchmarks

When adding performance-critical code:

```rust
use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn bench_my_function(c: &mut Criterion) {
    c.bench_function("my_function", |b| {
        b.iter(|| black_box(my_function()));
    });
}

criterion_group!(benches, bench_my_function);
criterion_main!(benches);
```

## 🐛 Reporting Bugs

When reporting bugs, please include:

- **Description:** Clear description of the bug
- **Steps to Reproduce:** Minimal steps to reproduce the issue
- **Expected Behavior:** What you expected to happen
- **Actual Behavior:** What actually happened
- **Environment:**
  - OS and version
  - Rust version (`rustc --version`)
  - Centrevo version
- **Additional Context:** Logs, screenshots, etc.

**Bug Report Template:**
```markdown
## Description
Brief description of the bug

## Steps to Reproduce
1. Step one
2. Step two
3. Step three

## Expected Behavior
What should happen

## Actual Behavior
What actually happens

## Environment
- OS: [e.g., Ubuntu 22.04]
- Rust: [e.g., 1.75.0]
- Centrevo: [e.g., 0.1.0]

## Additional Context
Any other relevant information
```

## 💡 Suggesting Features

When suggesting features:

- **Check existing issues** to avoid duplicates
- **Describe the use case** and why it's valuable
- **Provide examples** of how it would be used
- **Consider backwards compatibility**

**Feature Request Template:**
```markdown
## Problem/Need
What problem does this solve?

## Proposed Solution
How should it work?

## Alternatives Considered
Other approaches you've thought about

## Additional Context
Examples, mockups, related issues, etc.
```

## 📚 Documentation

Good documentation is crucial:

- **Code comments:** Explain *why*, not just *what*
- **API docs:** Use rustdoc for all public APIs
- **User guides:** Update CLI.md, PYTHON.md, etc.
- **Examples:** Add to `examples/` directory
- **README:** Keep README.md up to date

## 🧪 Testing Guidelines

### Unit Tests

- Test individual functions and methods
- Use `#[cfg(test)]` modules
- Mock external dependencies when needed

### Integration Tests

- Test complete workflows
- Place in `tests/` directory
- Test CLI commands end-to-end

### Property-Based Tests

For complex logic, consider property-based tests:
```rust
use quickcheck::QuickCheck;

#[test]
fn test_property() {
    fn prop(input: Vec<u8>) -> bool {
        // Property that should hold
        true
    }
    QuickCheck::new().quickcheck(prop as fn(Vec<u8>) -> bool);
}
```

## 🔍 Code Review

All submissions require review. We expect:

- **Timely responses** to feedback
- **Professional, respectful** communication
- **Willingness to iterate** on solutions
- **Learning from feedback**

As a reviewer:

- Be constructive and specific
- Explain *why* changes are needed
- Praise good solutions
- Be patient with new contributors

## 🏷️ Commit Messages

Good commit messages help maintain project history:

**Format:**
```
<type>(<scope>): <subject>

<body>

<footer>
```

**Types:**
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `style`: Code style changes (formatting, etc.)
- `refactor`: Code refactoring
- `perf`: Performance improvements
- `test`: Adding or updating tests
- `chore`: Maintenance tasks

**Examples:**
```
feat(simulation): add variable mutation rate support

Implement position-specific mutation rates to allow modeling
of hypermutable regions in centromeric arrays.

Closes #123
```

```
fix(storage): handle database locking correctly

Fix issue where concurrent queries could cause deadlock
by properly setting WAL mode and busy timeout.

Fixes #456
```

## 📊 Performance Considerations

When contributing performance-critical code:

1. **Benchmark before and after:**
   ```bash
   cargo bench
   ```

2. **Profile if needed:**
   ```bash
   cargo flamegraph --bench my_benchmark
   ```

3. **Document performance characteristics:**
   - Time complexity
   - Space complexity
   - Scaling behavior

4. **Consider:**
   - Memory allocations
   - Cache locality
   - Parallelization opportunities

## 🎯 Good First Issues

New contributors should look for issues labeled:
- `good first issue`: Suitable for beginners
- `help wanted`: Extra attention needed
- `documentation`: Documentation improvements

## 🙋 Getting Help

- **Questions:** Open a [GitHub Discussion](https://github.com/YOUR_USERNAME/centrevo/discussions)
- **Chat:** Join our [Discord/Slack/Matrix] (if applicable)
- **Email:** Contact maintainers directly for sensitive issues

## 📜 Code of Conduct

This project adheres to a Code of Conduct. By participating, you are expected to uphold this code:

- **Be respectful** and inclusive
- **Be collaborative** and constructive
- **Be patient** with new contributors
- **Focus on what is best** for the community

## ✅ Checklist Before Submitting

- [ ] Code follows style guidelines
- [ ] All tests pass (`cargo test`)
- [ ] No clippy warnings (`cargo clippy`)
- [ ] Code is formatted (`cargo fmt`)
- [ ] Documentation is updated
- [ ] Tests are added for new functionality
- [ ] Commit messages are clear
- [ ] Branch is up to date with main

## 🎉 Recognition

Contributors are recognized in:
- README.md acknowledgments
- Git history
- Release notes
- (Future) Contributors page

Thank you for contributing to Centrevo! 🚀

---

**Questions?** Feel free to ask in [GitHub Discussions](https://github.com/YOUR_USERNAME/centrevo/discussions) or open an issue.
