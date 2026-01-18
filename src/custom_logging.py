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
        #TODO: consider adding log_level here?
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
    
    def section(self, title: str, log_level: str = 'INFO'):
        """Print section header with indentation.
        
        Args:
            title: Section title
        """
        indent = '  ' * self.indent_level
        separator = '=' * 60

        # Check if the level would be logged
        numeric_level = getattr(logging, log_level.upper(), logging.INFO)
        if self.logger.isEnabledFor(numeric_level):
            print(f"\n{indent}{separator}")
            print(f"{indent}{title}")
            print(f"{indent}{separator}")
        
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


def init_logging(log_level: str, config_file: str = 'logging_config.yaml') -> None:
    """Initialize logging system with specified level.
    
    Sets up logging infrastructure from YAML config and configures the
    root logger to use the specified level. Should be called once at
    application startup before any logging occurs.
    
    Args:
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        config_file: Path to logging configuration YAML file (default: logging_config.yaml)
        
    Raises:
        ValueError: If log_level is not a valid logging level
        
    Example:
        >>> init_logging('DEBUG')
        >>> logger = get_logger(__name__)
        >>> logger.debug("Debug message")  # Will be shown
        
        >>> init_logging('INFO')
        >>> logger.debug("Debug message")  # Will NOT be shown
    """
    # Validate log level
    #TODO: move to input file validator?
    valid_levels = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
    log_level_upper = log_level.upper()
    
    if log_level_upper not in valid_levels:
        raise ValueError(
            f"Invalid log level: '{log_level}'. "
            f"Must be one of: {', '.join(valid_levels)}"
        )
    
    # Load logging infrastructure from YAML
    config_path = Path(config_file)
    
    if config_path.exists():
        try:
            with open(config_path, 'r') as f:
                config = yaml.safe_load(f)
                logging.config.dictConfig(config['logging'])
            
        except Exception as e:
            # Fall back to basic config if YAML is malformed
            print(f"⚠️  Error loading {config_file}: {e}")
            print(f"⚠️  Using basic logging configuration")
            logging.basicConfig(
                level=logging.INFO,
                format='%(message)s',
                force=True
            )
    else:
        # YAML not found - use basic config
        print(f"⚠️  {config_file} not found, using basic configuration")
        logging.basicConfig(
            level=logging.INFO,
            format='%(message)s',
            force=True
        )
    
    # Set the root logger level (affects all loggers)
    numeric_level = getattr(logging, log_level_upper)
    logging.getLogger().setLevel(numeric_level)


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