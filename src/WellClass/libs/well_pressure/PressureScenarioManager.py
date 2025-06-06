from dataclasses import dataclass, field
from typing import Dict
from .PressureScenario import PressureScenario
import pandas as pd



@dataclass
class PressureScenarioManager:
    scenarios: Dict[str, PressureScenario] = field(default_factory=dict)

    def create_scenario(self, name: str, **kwargs) -> PressureScenario:
        # Create a new PressureScenario instance with the given name and parameters
        if 'init_curves' in kwargs:
            kwargs['init_curves'] = kwargs['init_curves'].copy(deep=True)

        if name in self.scenarios:
            raise ValueError(f"Scenario with name '{name}' already exists.")
        
        scenario = PressureScenario(name=name, **kwargs)
        

        # Additional logic to initialize the scenario based on kwargs
        self.scenarios[name] = scenario
        return scenario

    def compute_all_scenarios(self):
        # Iterate over all scenarios and compute their pressure profiles
        for scenario in self.scenarios.values():
            scenario.compute_pressure_profile()



    def get_scenarios_summary(self) -> pd.DataFrame:
        # Collect parameters for each scenario into a list of dictionaries
        scenarios_data = []
        for name, scenario in self.scenarios.items():
            scenario_data = {
                'name': name,
                'from_resrvr': scenario.from_resrvr,
                'z_MSAD': scenario.z_MSAD,
                'p_MSAD': scenario.p_MSAD,
                'z_MSAD_brine': scenario.z_MSAD_brine,
                'p_MSAD_brine': scenario.p_MSAD_brine,
                'z_resrv': scenario.z_resrv,
                'p_resrv': scenario.p_resrv,
                'z_fluid_contact': scenario.z_fluid_contact,
                'p_fluid_contact': scenario.p_fluid_contact,
                'p_delta': scenario.p_delta,
                'fluid_type': scenario.fluid_type,
                'fluid_composition': scenario.fluid_composition,
                'specific_gravity': scenario.specific_gravity,
            }
            scenarios_data.append(scenario_data)
        
        # Convert the list of dictionaries into a DataFrame
        scenarios_df = pd.DataFrame(scenarios_data)
        return scenarios_df

    def collate_scenario_profiles(self, common_data: pd.DataFrame) -> pd.DataFrame:
        scenario_data = []

        # Collate scenario-specific profiles
        for name, scenario in self.scenarios.items():
            profiles_df = pd.DataFrame({
                (name, 'fluid_pressure'): scenario.init_curves.fluid_pressure,
                (name, 'brine_pressure'): scenario.init_curves.brine_pressure,
            }, index=common_data.index)
            scenario_data.append(profiles_df)

        # Concatenate all scenario DataFrames along the column axis with MultiIndex columns
        scenario_profiles = pd.concat(scenario_data, axis=1)
        # Add common data to the MultiIndex DataFrame
        scenario_profiles = pd.concat([common_data, scenario_profiles], axis=1)
        
        return scenario_profiles