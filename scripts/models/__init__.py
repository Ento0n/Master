"""Model definitions used by the interaction training pipeline."""

from .dscript import DScriptInteractionModel
from .fully_connected import FullyConnectedInteractionModel

__all__ = ["DScriptInteractionModel", "FullyConnectedInteractionModel"]
