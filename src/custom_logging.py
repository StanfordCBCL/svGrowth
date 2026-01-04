"""
Logging utilities for svGrowth G&R simulations.

Provides enhanced logging with indentation support and clean integration
with Python's standard logging module.
"""

import logging
import logging.config
import yaml
from pathlib import Path
from typing import Optional


class IndentedLogger:
    """Wrapper around Python logging with indentation support.
    
    Provides:
    - Indentation support (indent/dedent methods)
    - Context manager for automatic indent/dedent
    - Emoji prefixes for errors/warnings
    - Section headers
    
    Usage:
        logger = get_logger(__name__)
        
        logger.info("Starting computation")
        with logger:  # Auto-indent
            logger.debug("Step 1")
            logger.debug("Step 2")
        # Auto-dedent
        logger.info("Complete")
    """
    
    def __init__(self, name: str):
        """Initialize indented logger.
        
        Args:
            name: Logger name (typically __name__ from calling module)
        """
        self.logger = logging.getLogger(name)
        self.indent_level = 0
    
    def _format_msg(self, msg: str) -> str:
        """Add indentation to message.
        
        Args:
            msg: Message to format
            
        Returns:
            Indented message string
        """
        indent = '  ' * self.indent_level
        return f"{indent}{msg}"
    
    # ========================================================================
    # Standard logging levels
    # ========================================================================
    
    def debug(self, msg: str):
        """Log debug message with indentation.
        
        Args:
            msg: Debug message
        """
        self.logger.debug(self._format_msg(msg))
    
    def info(self, msg: str):
        """Log info message with indentation.
        
        Args:
            msg: Info message
        """
        self.logger.info(self._format_msg(msg))
    
    def warning(self, msg: str):
        """Log warning with emoji and indentation.
        
        Args:
            msg: Warning message
        """
        self.logger.warning(self._format_msg(f"⚠️  {msg}"))
    
    def error(self, msg: str):
        """Log error with emoji and indentation.
        
        Args:
            msg: Error message
        """
        self.logger.error(self._format_msg(f"❌ {msg}"))
    
    def critical(self, msg: str):
        """Log critical error with emoji and indentation.
        
        Args:
            msg: Critical error message
        """
        self.logger.critical(self._format_msg(f"🔥 {msg}"))
    
    # ========================================================================
    # Special formatting
    # ========================================================================
    
    def section(self, title: str):
        """Print section header with indentation.
        
        Args:
            title: Section title
        """
        indent = '  ' * self.indent_level
        separator = '=' * 60
        self.logger.info(f"\n{indent}{separator}")
        self.logger.info(f"{indent}{title}")
        self.logger.info(f"{indent}{separator}")
    
    # ========================================================================
    # Indentation control
    # ========================================================================
    
    def indent(self):
        """Increase indentation level."""
        self.indent_level += 1
    
    def dedent(self):
        """Decrease indentation level."""
        self.indent_level = max(0, self.indent_level - 1)
    
    # ========================================================================
    # Context manager support
    # ========================================================================
    
    def __enter__(self):
        """Enter indented context (auto-indent).
        
        Returns:
            Self for context manager
        """
        self.indent()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Exit indented context (auto-dedent).
        
        Args:
            exc_type: Exception type (if raised)
            exc_val: Exception value (if raised)
            exc_tb: Exception traceback (if raised)
            
        Returns:
            False to propagate exceptions
        """
        self.dedent()
        return False  # Don't suppress exceptions


def setup_logging(config_file: str = 'logging_config.yaml') -> None:
    """Initialize logging system from YAML configuration file.
    
    Should be called once at application startup before any logging occurs.
    Falls back to basic configuration if file is not found or invalid.
    
    Args:
        config_file: Path to logging configuration YAML file
        
    Example:
        >>> setup_logging('logging_config.yaml')
        >>> logger = get_logger(__name__)
        >>> logger.info("Application started")
    """
    config_path = Path(config_file)
    
    # Try to load from YAML config
    if config_path.exists():
        try:
            with open(config_path, 'r') as f:
                config = yaml.safe_load(f)
                logging.config.dictConfig(config['logging'])
            
            # Log successful initialization
            root_logger = logging.getLogger()
            root_logger.info(f"✓ Logging initialized from {config_file}")
            return
            
        except Exception as e:
            # Fall through to basic config if YAML is malformed
            print(f"⚠️  Error loading {config_file}: {e}")
    else:
        print(f"⚠️  {config_file} not found")
    
    # Fallback: basic configuration
    logging.basicConfig(
        level=logging.INFO,
        format='%(message)s',
        force=True  # Override any existing config
    )
    
    root_logger = logging.getLogger()
    root_logger.warning("Using basic logging configuration")


def get_logger(name: str) -> IndentedLogger:
    """Get enhanced logger with indentation support.
    
    Args:
        name: Logger name (typically __name__ from calling module)
        
    Returns:
        IndentedLogger instance
        
    Example:
        >>> logger = get_logger(__name__)
        >>> logger.info("Starting process")
        >>> with logger:
        ...     logger.debug("Detailed step")
        >>> logger.info("Process complete")
    """
    return IndentedLogger(name)