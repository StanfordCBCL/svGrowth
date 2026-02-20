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


class Box:
    """Handles box-drawing formatting for structured output.
    
    Box structure visualization:
        ┌─ Title ─────────────────────  ← box_start()
        │  Content line                  ← format_line(spacing=2)
        │    Indented content            ← format_line(spacing=4)
        │    Label:      Value           ← format_line(spacing=4, aligned)
        └─────────────────────────────  ← box_end()
    """
    
    # Box drawing characters
    VERTICAL = "│"
    HORIZONTAL = "─"
    TOP_LEFT = "┌"
    BOTTOM_LEFT = "└"
    
    def __init__(self, indent: str = "", outer_box_prefix: str = ""):
        """Initialize box formatter.
        
        Args:
            indent: Base indentation string (e.g., "" for layer, "  " for constituent)
            outer_box_prefix: Prefix from outer boxes (e.g., "│  " for nested boxes)
        """
        self.indent = indent
        self.outer_box_prefix = outer_box_prefix
    
    def box_start(self, title: str, width: int = 60) -> str:
        """Format box opening with title.
        
        Args:
            title: Box title text
            width: Total width of the box border
            
        Returns:
            Formatted box header string
            
        Visual output:
            ┌─ Title ─────────────────────────────────────
        """
        # Calculate separator to match box_end width exactly
        # Total: outer_prefix + indent + "┌─ " + title + " " + separator
        # We want total line length = outer_prefix + indent + width
        used_chars = len(f"┌─ {title} ")
        separator_length = width - used_chars
        separator = self.HORIZONTAL * separator_length
        
        return f"{self.outer_box_prefix}{self.indent}{self.TOP_LEFT}{self.HORIZONTAL} {title} {separator}"
    
    def box_end(self, width: int = 60) -> str:
        """Format box closing.
        
        Args:
            width: Total width of the box border
            
        Returns:
            Formatted box footer string
            
        Visual output:
            └─────────────────────────────────────────────
        """
        # Match box_start width exactly
        separator = self.HORIZONTAL * (width - 1)  # -1 for the └ character
        return f"{self.outer_box_prefix}{self.indent}{self.BOTTOM_LEFT}{separator}"
    
    def format_line(self, text: str = "", spacing: int = 2, label: str = None, label_width: int = 12) -> str:
        """Format a line within the box.
        
        Args:
            text: Content to display (or value if label is provided)
            spacing: Number of spaces after vertical bar (default: 2)
            label: Optional label for aligned key-value pairs
            label_width: Width for label column when using aligned mode
            
        Returns:
            Formatted line string
            
        Visual output examples:
            format_line("Content", spacing=2)
            → │  Content
            
            format_line("Item", spacing=4)
            → │    Item
            
            format_line("Value", spacing=4, label="Key:", label_width=12)
            → │    Key:         Value
            
            format_line("", spacing=0)
            → │
        """
        if not text and label is None:
            # Blank line
            return f"{self.outer_box_prefix}{self.indent}{self.VERTICAL}"
        
        spaces = " " * spacing
        
        if label is not None:
            # Aligned key-value pair
            return f"{self.outer_box_prefix}{self.indent}{self.VERTICAL}{spaces}{label:<{label_width}} {text}"
        else:
            # Plain content
            return f"{self.outer_box_prefix}{self.indent}{self.VERTICAL}{spaces}{text}"


class BoxContext:
    """Context manager for box-drawing mode.
    
    Handles box lifecycle: creates box, prints header, manages indentation,
    prints footer, and restores previous box state.
    
    Usage:
        with logger.box("Layer: Cerebral artery"):
            logger.box_line("Kinematics: ThinWall")
    """
    
    def __init__(self, logger: 'CustomLogger', title: str, width: int = 60):
        """Initialize box context.
        
        Args:
            logger: Parent CustomLogger instance
            title: Title for box header
            width: Width of box border
        """
        self.logger = logger
        self.title = title
        self.width = width
        self.previous_box = None
    
    def __enter__(self):
        """Enter box context - print header and create box.
        
        Creates new box at current indentation level, prints header,
        saves previous box, and increases indentation for box contents.
        """
        # Create new box with current indentation
        indent = '  ' * self.logger.indent_level
        
        # Build outer box prefix (continuation of parent boxes)
        outer_prefix = ""
        if self.logger.indent_level > 0:
            # We're inside another box - need to continue its vertical line
            indent = "" # Don't add extra indent for nested box - it will be handled by outer box prefix
            parent_box = self.logger.current_box
            outer_prefix = f"{parent_box.outer_box_prefix}{parent_box.indent}{Box.VERTICAL}  "
        
        # Create new box with current indentation and outer prefix
        new_box = Box(indent=indent, outer_box_prefix=outer_prefix)
        
        # Print box opening
        header = new_box.box_start(self.title, self.width)
        self.logger.logger.info(header)
        
        # Save previous box and set new one
        self.previous_box = self.logger.current_box
        self.logger.current_box = new_box
        
        # Increase indent for content inside box
        self.logger.indent()
        
        return self.logger
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Exit box context - print footer and restore previous box.
        
        Decreases indentation, prints box footer, and restores the
        previous box state.
        """
        # Decrease indent
        self.logger.dedent()
        
        # Print box closing
        footer = self.logger.current_box.box_end(self.width)
        self.logger.logger.info(footer)
        
        # Restore previous box
        self.logger.current_box = self.previous_box
        
        return False


class CustomLogger:
    """Wrapper around Python logging with indentation and box-drawing support.
    
    Provides:
    - Indentation support (indent/dedent methods)
    - Context manager for automatic indent/dedent
    - Box-drawing utilities for structured output
    - Emoji prefixes for errors/warnings
    - Section headers
    
    Usage:
        logger = get_logger(__name__)
        
        # Standard logging with auto-indentation
        logger.info("Starting computation")
        with logger:
            logger.debug("Step 1")
            logger.debug("Step 2")
        logger.info("Complete")
        
        # Box-drawing for structured output
        with logger.box("Layer: Arterial Wall"):
            logger.box_line("Kinematics: ThinWall")
            logger.box_section("Geometry")
            logger.box_item("Inner radius: 1.4 mm")
    """
    
    def __init__(self, name: str):
        """Initialize indented logger.
        
        Args:
            name: Logger name (typically __name__ from calling module)
        """
        #TODO: consider adding log_level here?
        self.logger = logging.getLogger(name)
        self.indent_level = 0
        self.current_box = Box(indent="")  # Default box for when not in box context
    
    # =========================================================================
    # BOX CONTEXT MANAGER
    # =========================================================================
    
    def box(self, title: str, width: int = 60):
        """Create box context manager for structured output.
        
        Args:
            title: Title for box header
            width: Width of box border (default: 60)
            
        Returns:
            BoxContext for use with 'with' statement
            
        Example:
            with logger.box("Layer: Cerebral artery"):
                logger.box_line("Kinematics: ThinWall")
                logger.box_section("Homeostatic Geometry")
                logger.box_item("Inner radius: 1.4 mm")
        """
        return BoxContext(self, title, width)
    
    # =========================================================================
    # BOX DRAWING METHODS
    # =========================================================================
    
    def box_line(self, text: str = ""):
        """Print content line within current box.
        
        Args:
            text: Content to display (empty for blank line)
            
        Visual output:
            │  Text content here
            │  (blank if text is empty)
        """
        line = self.current_box.format_line(text, spacing=2)
        self.logger.info(line)
    
    def box_section(self, title: str):
        """Print section header within current box.
        
        Adds blank line before section title for visual separation.
        
        Args:
            title: Section title
            
        Visual output:
            │
            │  Section Title
        """
        self.logger.info(self.current_box.format_line("", spacing=0))  # Blank line
        self.logger.info(self.current_box.format_line(title, spacing=2))  # Section title
    
    def box_item(self, text: str, extra_indent: int = 2):
        """Print indented item within a section.
        
        Args:
            text: Item content
            extra_indent: Additional indentation spaces (default: 2)
            
        Visual output (extra_indent=2):
            │    Item content here
        """
        item = self.current_box.format_line(text, spacing=2 + extra_indent)
        self.logger.info(item)
    
    def box_item_aligned(self, label: str, value: str, label_width: int = 12, extra_indent: int = 2):
        """Print aligned key-value pair within box.
        
        Args:
            label: Left-aligned label
            value: Value to display
            label_width: Width for label column (default: 12)
            extra_indent: Additional indentation spaces (default: 2)
            
        Visual output (label_width=12, extra_indent=2):
            │    Label:       Value here
        """
        item = self.current_box.format_line(
            text=value, 
            spacing=2 + extra_indent, 
            label=label, 
            label_width=label_width
        )
        self.logger.info(item)
    
    # =========================================================================
    # STANDARD LOGGING METHODS
    # =========================================================================
    
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
    
    # =========================================================================
    # FORMATTING METHODS
    # =========================================================================
    
    def section(self, title: str, log_level: str = 'INFO'):
        """Print section header with separator lines.
        
        Args:
            title: Section title
            log_level: Logging level (default: 'INFO')
            
        Visual output:
            ============================================================
            Section Title
            ============================================================
        """
        indent = '  ' * self.indent_level
        separator = '=' * 60

        # Check if the level would be logged
        numeric_level = getattr(logging, log_level.upper(), logging.INFO)
        if self.logger.isEnabledFor(numeric_level):
            print(f"\n{indent}{separator}")
            print(f"{indent}{title}")
            print(f"{indent}{separator}")
    
    # =========================================================================
    # INDENTATION CONTROL
    # =========================================================================
    
    def indent(self):
        """Increase indentation level by one step (2 spaces)."""
        self.indent_level += 1
    
    def dedent(self):
        """Decrease indentation level by one step (2 spaces)."""
        self.indent_level = max(0, self.indent_level - 1)

    def _format_msg(self, msg: str) -> str:
        """Add indentation to message.
        
        Args:
            msg: Message to format
            
        Returns:
            Indented message string
        """
        indent = '  ' * self.indent_level
        return f"{indent}{msg}"
    
    # =========================================================================
    # CONTEXT MANAGER (for simple indentation)
    # =========================================================================
    
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
    # TODO: config file currently don't work if not in root dir. Consider making path handling more robust.
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
    numeric_level = getattr(logging, log_level.upper())
    logging.getLogger().setLevel(numeric_level)


def get_logger(name: str) -> CustomLogger:
    """Get enhanced logger with indentation and box-drawing support.
    
    Args:
        name: Logger name (typically __name__ from calling module)
        
    Returns:
        CustomLogger instance
        
    Example:
        >>> logger = get_logger(__name__)
        >>> logger.info("Starting process")
        >>> with logger:
        ...     logger.debug("Detailed step")
        >>> logger.info("Process complete")
        
        >>> with logger.box("Configuration"):
        ...     logger.box_line("Parameter: value")
    """
    return CustomLogger(name)