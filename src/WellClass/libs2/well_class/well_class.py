import json
from dataclasses import dataclass
from pathlib import Path

from src.WellClass.libs2.models.well_model import WellModel
from src.WellClass.libs2.well_class.well_raw import WellRaw


@dataclass
class Well(WellRaw):
    def __post_init__(self) -> None:
        """Compute basic well information"""
        self._check_inventory()

    def _check_inventory(self) -> None:
        self.inventory = {
            "hole_casings": bool(self.hole_casings),
            "survey": bool(self.survey),
            "plugs": bool(self.plugs),
            "stratigraphy": bool(self.stratigraphy),
        }

    @classmethod
    def from_pydantic(cls, model: WellModel) -> "Well":
        return cls(
            header=model.spec.well_header.model_dump(),
            hole_casings=[hc.model_dump() for hc in model.spec.hole_casings] if model.spec.hole_casings else None,
            survey=model.spec.well_survey.model_dump() if model.spec.well_survey else None,
            plugs=[pl.model_dump() for pl in model.spec.plugs] if model.spec.plugs else None,
            stratigraphy=[st.model_dump() for st in model.spec.stratigraphy] if model.spec.stratigraphy else None,
        )

    @classmethod
    def from_json(cls, json_file: str | Path) -> "Well":
        if isinstance(json_file, str):
            json_file = Path(json_file)
        with json_file.open(encoding="utf-8") as f:
            json_data = json.load(f)

        model = WellModel.model_validate(json_data)
        return cls.from_pydantic(model)
