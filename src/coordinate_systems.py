"""
Coordinate system utilities for svGrowth.

Centralizes all coordinate system information with type-safe enums:
- Direction names (radial, circumferential, axial)
- Tensor indices ([0,0], [1,1], [2,2] and off-diagonal)
- Fiber angle mappings
- Variable naming conventions

Makes code coordinate-system agnostic and extensible.
"""

from enum import Enum
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import numpy as np


# =============================================================================
# COORDINATE SYSTEM TYPE ENUM
# =============================================================================

class CoordinateSystemType(Enum):
    """Supported coordinate systems."""
    CYLINDRICAL = 'cylindrical'
    CARTESIAN = 'cartesian'
    SPHERICAL = 'spherical'


# =============================================================================
# DIRECTION ENUMS (Type-safe for each coordinate system)
# =============================================================================

class CylindricalDirection(Enum):
    """Directions in cylindrical coordinates."""
    RADIAL = 'radial'
    CIRCUMFERENTIAL = 'circumferential'
    AXIAL = 'axial'


class CartesianDirection(Enum):
    """Directions in Cartesian coordinates."""
    X = 'x-direction'
    Y = 'y-direction'
    Z = 'z-direction'


class SphericalDirection(Enum):
    """Directions in spherical coordinates."""
    RADIAL = 'radial'
    POLAR = 'polar'
    AZIMUTHAL = 'azimuthal'


# =============================================================================
# DIRECTION INFO
# =============================================================================

@dataclass(frozen=True)
class DirectionInfo:
    """Information about a coordinate direction.
    
    Represents a basis direction in a coordinate system (e.g., radial, circumferential).
    NOT a tensor component (which would be an (i,j) pair in a tensor).
    
    Examples:
        Cylindrical radial:
            name='radial', index=0, symbol='r', variable_subscript='r'
        
        Cylindrical circumferential:
            name='circumferential', index=1, symbol='θ', variable_subscript='theta'
        
        Cylindrical axial:
            name='axial', index=2, symbol='z', variable_subscript='z'
    """
    name: str                   # Human-readable name ('radial', 'circumferential', 'axial')
    index: int                  # Position in coordinate system (0, 1, or 2)
    symbol: str                 # Math symbol ('r', 'θ', 'z')
    variable_subscript: str     # For variable naming ('r', 'theta', 'z')
    
    def get_diagonal_index(self) -> Tuple[int, int]:
        """Get diagonal tensor index for this direction."""
        return (self.index, self.index)
    
    def get_variable_name(self, quantity: str) -> str:
        """Get variable name for this direction.
        
        Args:
            quantity: Base quantity name ('lambda', 'sigma', 'G', etc.)
            
        Returns:
            Full variable name
            
        Examples:
            >>> dir_info.get_variable_name('lambda')
            'lambda_theta'
            >>> dir_info.get_variable_name('sigma')
            'sigma_theta'
            >>> dir_info.get_variable_name('G')
            'G_theta'
        """
        return f"{quantity}_{self.variable_subscript}"
    
    # Convenience methods for common quantities
    def get_stretch_name(self) -> str:
        """Get stretch variable name: lambda_r, lambda_theta, etc."""
        return self.get_variable_name('lambda')
    
    def get_stress_name(self) -> str:
        """Get stress variable name: sigma_r, sigma_theta, etc."""
        return self.get_variable_name('sigma')
    
    def get_deposition_stretch_name(self) -> str:
        """Get deposition stretch variable name: G_r, G_theta, etc."""
        return self.get_variable_name('G')


# =============================================================================
# COORDINATE SYSTEM METADATA (Registry)
# =============================================================================

COORDINATE_SYSTEM_INFO = {
    'cylindrical': {
        'description': 'Cylindrical coordinates (r, θ, z) for thin-wall vessels',
        'directions': ['radial', 'circumferential', 'axial'],
        'symbols': ['r', 'θ', 'z'],
    },
    'cartesian': {
        'description': 'Cartesian coordinates (x, y, z) for general 3D problems',
        'directions': ['x-direction', 'y-direction', 'z-direction'],
        'symbols': ['x', 'y', 'z'],
    },
    'spherical': {
        'description': 'Spherical coordinates (r, θ, φ) for spherical geometries',
        'directions': ['radial', 'polar', 'azimuthal'],
        'symbols': ['r', 'θ', 'φ'],
    }
}


# =============================================================================
# COORDINATE SYSTEM CLASS
# =============================================================================

class CoordinateSystem:
    """A coordinate system with direction information and utilities.
    
    Provides type-safe mapping between:
    - Physical directions (CylindricalDirection.RADIAL, etc.)
    - Tensor indices ([0,0], [1,1], [2,2] and off-diagonal)
    - Variable names (sigma_theta, G_r, lambda_z)
    - Fiber angles (0° = axial, 90° = circumferential)
    
    Usage:
        >>> coord_system = CoordinateSystem.get(CoordinateSystemType.CYLINDRICAL)
        
        >>> # Type-safe direction access
        >>> theta_dir = coord_system.get_direction(CylindricalDirection.CIRCUMFERENTIAL)
        >>> print(theta_dir.index)  # 1
        >>> print(theta_dir.symbol)  # 'θ'
        
        >>> # Build diagonal tensor
        >>> stretches = {
        ...     CylindricalDirection.RADIAL: 0.7,
        ...     CylindricalDirection.CIRCUMFERENTIAL: 1.4,
        ...     CylindricalDirection.AXIAL: 1.0
        ... }
        >>> F = coord_system.build_diagonal_tensor(stretches)
        
        >>> # Extract value
        >>> lambda_theta = coord_system.get_diagonal_value(
        ...     F, CylindricalDirection.CIRCUMFERENTIAL
        ... )
    """
    
    # Singleton cache
    _instances: Dict[CoordinateSystemType, 'CoordinateSystem'] = {}
    
    def __init__(self, system_type: CoordinateSystemType):
        """Initialize coordinate system.
        
        Args:
            system_type: Type of coordinate system
        """
        self.system_type = system_type
        self.directions: Dict[Enum, DirectionInfo] = {}  # Enum → DirectionInfo
        self._direction_list: List[DirectionInfo] = []    # Ordered list for indexing
        self._fiber_angle_to_direction: Dict[float, int] = {}
        
        # Initialize based on type
        if system_type == CoordinateSystemType.CYLINDRICAL:
            self.direction_enum = CylindricalDirection
            self._init_cylindrical()
        elif system_type == CoordinateSystemType.CARTESIAN:
            self.direction_enum = CartesianDirection
            self._init_cartesian()
        elif system_type == CoordinateSystemType.SPHERICAL:
            self.direction_enum = SphericalDirection
            self._init_spherical()
        else:
            raise ValueError(f"Unknown coordinate system type: {system_type}")
    
    def _init_cylindrical(self):
        """Initialize cylindrical coordinates (r, θ, z)."""
        # Build directions with enum keys
        self.directions = {
            CylindricalDirection.RADIAL: 
                DirectionInfo(
                    name='radial',
                    index=0,
                    symbol='r',
                    variable_subscript='r'
                ),
            CylindricalDirection.CIRCUMFERENTIAL: 
                DirectionInfo(
                    name='circumferential',
                    index=1,
                    symbol='θ',
                    variable_subscript='theta'
                ),
            CylindricalDirection.AXIAL: 
                DirectionInfo(
                    name='axial',
                    index=2,
                    symbol='z',
                    variable_subscript='z'
                )
        }
        
        # Ordered list for index-based access
        self._direction_list = [
            self.directions[CylindricalDirection.RADIAL],
            self.directions[CylindricalDirection.CIRCUMFERENTIAL],
            self.directions[CylindricalDirection.AXIAL]
        ]
        
        # Fiber angle → direction mapping
        # Convention: angle measured from axial (z) direction
        self._fiber_angle_to_direction = {
            0.0: 2,    # 0° = axial (z-direction, index 2)
            90.0: 1,   # 90° = circumferential (θ-direction, index 1)
            # Future: 45.0 for helical (requires special handling)
        }
    
    def _init_cartesian(self):
        """Initialize Cartesian coordinates (x, y, z)."""
        self.directions = {
            CartesianDirection.X: 
                DirectionInfo(
                    name='x-direction',
                    index=0,
                    symbol='x',
                    variable_subscript='x'
                ),
            CartesianDirection.Y: 
                DirectionInfo(
                    name='y-direction',
                    index=1,
                    symbol='y',
                    variable_subscript='y'
                ),
            CartesianDirection.Z: 
                DirectionInfo(
                    name='z-direction',
                    index=2,
                    symbol='z',
                    variable_subscript='z'
                )
        }
        
        self._direction_list = [
            self.directions[CartesianDirection.X],
            self.directions[CartesianDirection.Y],
            self.directions[CartesianDirection.Z]
        ]
        
        # No simple fiber angle mapping for Cartesian (use direction vectors)
        self._fiber_angle_to_direction = {}
    
    def _init_spherical(self):
        """Initialize spherical coordinates (r, θ, φ)."""
        self.directions = {
            SphericalDirection.RADIAL: 
                DirectionInfo(
                    name='radial',
                    index=0,
                    symbol='r',
                    variable_subscript='r'
                ),
            SphericalDirection.POLAR: 
                DirectionInfo(
                    name='polar',
                    index=1,
                    symbol='θ',
                    variable_subscript='theta'
                ),
            SphericalDirection.AZIMUTHAL: 
                DirectionInfo(
                    name='azimuthal',
                    index=2,
                    symbol='φ',
                    variable_subscript='phi'
                )
        }
        
        self._direction_list = [
            self.directions[SphericalDirection.RADIAL],
            self.directions[SphericalDirection.POLAR],
            self.directions[SphericalDirection.AZIMUTHAL]
        ]
        
        # Fiber angle mapping TBD for spherical
        self._fiber_angle_to_direction = {}
    
    # =========================================================================
    # DIRECTION ACCESS (Type-safe!)
    # =========================================================================
    
    def get_direction(self, direction: Enum) -> DirectionInfo:
        """Get direction info by enum (type-safe).
        
        Args:
            direction: Direction enum (e.g., CylindricalDirection.RADIAL)
            
        Returns:
            DirectionInfo
            
        Raises:
            ValueError: If direction not valid for this coordinate system
            
        Example:
            >>> dir_info = coord_system.get_direction(CylindricalDirection.CIRCUMFERENTIAL)
            >>> dir_info.index  # 1
        """
        if direction not in self.directions:
            raise ValueError(
                f"Invalid direction {direction} for {self.system_type.value} coordinates. "
                f"Valid: {list(self.directions.keys())}"
            )
        return self.directions[direction]
    
    def get_direction_by_index(self, index: int) -> DirectionInfo:
        """Get direction by tensor index (0, 1, or 2).
        
        Args:
            index: Tensor index
            
        Returns:
            DirectionInfo
            
        Raises:
            ValueError: If index not in [0, 1, 2]
        """
        if not (0 <= index <= 2):
            raise ValueError(f"Index must be 0, 1, or 2, got {index}")
        return self._direction_list[index]
    
    # =========================================================================
    # TENSOR CONSTRUCTION (What kinematics needs!)
    # =========================================================================
    
    def build_diagonal_tensor(
        self,
        directions: Dict[Enum, float]
    ) -> np.ndarray:
        """Build diagonal tensor from direction values (type-safe).
        
        Args:
            directions: Dict mapping direction enum → value
                Example: {
                    CylindricalDirection.RADIAL: 0.7,
                    CylindricalDirection.CIRCUMFERENTIAL: 1.4,
                    CylindricalDirection.AXIAL: 1.0
                }
                
        Returns:
            3×3 diagonal tensor in correct order
            
        Raises:
            ValueError: If any direction is missing
        """
        # Validate all directions provided
        for dir_enum in self.directions.keys():
            if dir_enum not in directions:
                raise ValueError(
                    f"Missing direction {dir_enum}. "
                    f"Required: {list(self.directions.keys())}"
                )
    
        # Build tensor in correct order using ordered list
        # Use the enum from self.directions dict, not dir_info.name
        diagonal_values = []
        for dir_info in self._direction_list:
            # Find the enum key that maps to this dir_info
            for enum_key, info in self.directions.items():
                if info is dir_info:
                    diagonal_values.append(directions[enum_key])
                    break
    
        return np.diag(diagonal_values)
    
    def get_diagonal_value(
        self,
        tensor: np.ndarray,
        direction: Enum
    ) -> float:
        """Extract diagonal value for direction (type-safe).
        
        Args:
            tensor: 3×3 tensor
            direction: Direction enum
            
        Returns:
            Diagonal value
            
        Example:
            >>> lambda_theta = coord_system.get_diagonal_value(
            ...     F, CylindricalDirection.CIRCUMFERENTIAL
            ... )
        """
        dir_info = self.get_direction(direction)
        return tensor[dir_info.index, dir_info.index]
    
    def set_diagonal_value(
        self,
        tensor: np.ndarray,
        direction: Enum,
        value: float
    ) -> None:
        """Set diagonal value for direction (in-place, type-safe).
        
        Args:
            tensor: 3×3 tensor to modify
            direction: Direction enum
            value: New value
            
        Example:
            >>> coord_system.set_diagonal_value(
            ...     sigma, CylindricalDirection.RADIAL, 0.0
            ... )
        """
        dir_info = self.get_direction(direction)
        tensor[dir_info.index, dir_info.index] = value
    
    def extract_diagonal_directions(
        self,
        tensor: np.ndarray
    ) -> Dict[Enum, float]:
        """Extract all diagonal values as dict (type-safe).
        
        Args:
            tensor: 3×3 tensor
            
        Returns:
            Dict mapping direction enum → diagonal value
            
        Example:
            >>> stretches = coord_system.extract_diagonal_directions(F)
            >>> stretches[CylindricalDirection.CIRCUMFERENTIAL]  # 1.4
        """
        return {
            dir_enum: tensor[dir_info.index, dir_info.index]
            for dir_enum, dir_info in self.directions.items()
        }
    
    # =========================================================================
    # TENSOR INDEX ACCESS (including off-diagonal)
    # =========================================================================
    
    def get_direction_pair_by_indices(
        self, 
        i: int, 
        j: int
    ) -> Tuple[DirectionInfo, DirectionInfo]:
        """Get direction pair for tensor element [i, j].
        
        Args:
            i: Row index (0, 1, or 2)
            j: Column index (0, 1, or 2)
            
        Returns:
            Tuple of (row_direction, col_direction)
            
        Example:
            >>> dir_pair = coord_system.get_direction_pair_by_indices(0, 1)
            >>> # (DirectionInfo(radial), DirectionInfo(circumferential))
        """
        if not (0 <= i <= 2 and 0 <= j <= 2):
            raise ValueError(f"Indices must be 0, 1, or 2, got ({i}, {j})")
        
        return (self._direction_list[i], self._direction_list[j])
    
    def get_diagonal_indices(self) -> List[Tuple[int, int]]:
        """Get list of diagonal indices: [(0,0), (1,1), (2,2)]."""
        return [(i, i) for i in range(3)]
    
    def get_off_diagonal_indices(self) -> List[Tuple[int, int]]:
        """Get list of off-diagonal indices (upper triangle): [(0,1), (0,2), (1,2)]."""
        indices = []
        for i in range(3):
            for j in range(i + 1, 3):
                indices.append((i, j))
        return indices
    
    def get_all_indices(self) -> List[Tuple[int, int]]:
        """Get all unique tensor indices (upper triangle): 
        [(0,0), (0,1), (0,2), (1,1), (1,2), (2,2)]."""
        indices = []
        for i in range(3):
            for j in range(i, 3):
                indices.append((i, j))
        return indices
    
    # =========================================================================
    # FIBER ANGLE MAPPING
    # =========================================================================
    
    def get_tensor_index_from_fiber_angle(self, angle_degrees: float) -> int:
        """Get tensor index for fiber at given angle.
        
        Args:
            angle_degrees: Fiber angle in degrees
            
        Returns:
            Tensor index (0, 1, or 2)
            
        Raises:
            NotImplementedError: If angle not supported
        """
        # Check exact match
        if angle_degrees in self._fiber_angle_to_direction:
            return self._fiber_angle_to_direction[angle_degrees]
        
        # Check approximate match (within 0.1°)
        for supported_angle, index in self._fiber_angle_to_direction.items():
            if abs(angle_degrees - supported_angle) < 0.1:
                return index
        
        # Not supported
        supported_angles = list(self._fiber_angle_to_direction.keys())
        raise NotImplementedError(
            f"Fiber angle {angle_degrees}° not supported in "
            f"{self.system_type.value} coordinates. "
            f"Supported angles: {supported_angles}"
        )
    
    def get_direction_from_fiber_angle(self, angle_degrees: float) -> Enum:
        """Get direction enum for fiber at given angle.
        
        Args:
            angle_degrees: Fiber angle in degrees
            
        Returns:
            Direction enum (e.g., CylindricalDirection.CIRCUMFERENTIAL)
            
        Raises:
            NotImplementedError: If angle not supported
        """
        index = self.get_tensor_index_from_fiber_angle(angle_degrees)
        dir_info = self.get_direction_by_index(index)
        return dir_info.name
    
    # =========================================================================
    # VARIABLE NAMING
    # =========================================================================
    
    def get_all_stretch_names(self) -> List[str]:
        """Get all stretch variable names: ['lambda_r', 'lambda_theta', 'lambda_z']."""
        return [dir_info.get_stretch_name() for dir_info in self._direction_list]
    
    def get_all_stress_names(self) -> List[str]:
        """Get all stress variable names: ['sigma_r', 'sigma_theta', 'sigma_z']."""
        return [dir_info.get_stress_name() for dir_info in self._direction_list]
    
    def get_all_deposition_stretch_names(self) -> List[str]:
        """Get all deposition stretch variable names: ['G_r', 'G_theta', 'G_z']."""
        return [dir_info.get_deposition_stretch_name() for dir_info in self._direction_list]
    
    # =========================================================================
    # FORMATTING UTILITIES
    # =========================================================================
    
    def format_tensor_component_name(self, i: int, j: int) -> str:
        """Get name for tensor component [i, j].
        
        Args:
            i: Row index
            j: Column index
            
        Returns:
            Formatted component name
            
        Examples:
            (0, 0) → 'rr' (diagonal)
            (0, 1) → 'rθ' (off-diagonal, shear)
            (1, 2) → 'θz' (off-diagonal, shear)
        """
        dir_i, dir_j = self.get_direction_pair_by_indices(i, j)
        
        if i == j:
            # Diagonal: use double subscript
            return f"{dir_i.symbol}{dir_i.symbol}"
        else:
            # Off-diagonal: use mixed subscript
            return f"{dir_i.symbol}{dir_j.symbol}"
    
    def format_tensor_diagonal(
        self, 
        tensor: np.ndarray, 
        precision: int = 4
    ) -> str:
        """Format tensor diagonal with direction names.
        
        Args:
            tensor: 3x3 tensor
            precision: Decimal precision
            
        Returns:
            Formatted string like "r: 1.0000, θ: 2.0000, z: 0.5000"
        """
        parts = []
        for dir_info in self._direction_list:
            value = tensor[dir_info.index, dir_info.index]
            parts.append(f"{dir_info.symbol}: {value:.{precision}f}")
        return ", ".join(parts)
    
    def format_full_tensor(
        self, 
        tensor: np.ndarray, 
        precision: int = 4,
        include_off_diagonal: bool = False
    ) -> str:
        """Format tensor with direction names.
        
        Args:
            tensor: 3x3 tensor
            precision: Decimal precision
            include_off_diagonal: Whether to include off-diagonal terms
            
        Returns:
            Formatted string
        """
        if include_off_diagonal:
            # Full tensor (3 rows)
            lines = []
            for i in range(3):
                row_parts = []
                for j in range(3):
                    name = self.format_tensor_component_name(i, j)
                    value = tensor[i, j]
                    row_parts.append(f"{name}: {value:.{precision}f}")
                lines.append(" | ".join(row_parts))
            return "\n".join(lines)
        else:
            # Diagonal only (current behavior)
            return self.format_tensor_diagonal(tensor, precision)
    
    # =========================================================================
    # FACTORY METHOD (Singleton)
    # =========================================================================
    
    @classmethod
    def get(cls, system_type: CoordinateSystemType) -> 'CoordinateSystem':
        """Get coordinate system instance (singleton pattern).
        
        Args:
            system_type: Type of coordinate system
            
        Returns:
            CoordinateSystem instance (reused if already created)
        """
        if system_type not in cls._instances:
            cls._instances[system_type] = cls(system_type)
        return cls._instances[system_type]
    
    # =========================================================================
    # STRING REPRESENTATION
    # =========================================================================
    
    def __str__(self) -> str:
        """String representation."""
        dir_names = [dir_info.name for dir_info in self._direction_list]
        return f"{self.system_type.value.capitalize()} ({', '.join(dir_names)})"
    
    def __repr__(self) -> str:
        """Developer representation."""
        return f"CoordinateSystem(type={self.system_type.value})"