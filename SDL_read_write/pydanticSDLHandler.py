################################################################################
### mtrk project - SDL file generator allowing to simply define MRI pulse    ###
### sequences that can be read by the mtrk project simulator and driver      ###
### sequence.                                                                ###
### Version 0.0.1                                                            ###
### Anais Artiges and the mtrk project team at NYU - 09/07/2023              ###
################################################################################  

from __future__ import annotations

from typing import List, Optional, Union
from typing_extensions import Literal, Any 
from pydantic import BaseModel, SerializeAsAny, Extra, root_validator
from decimal import Decimal

### file section
class File(BaseModel):
    format: str = "mtrk-SDL"
    version: int = 9999
    measurement: str = "default_measurement"
    system: str = "default_system"


### settings section
class Settings(BaseModel):
    readout_os: int = 9999

    class Config:
        extra = Extra.allow

    @root_validator(pre=True)
    def check_extra_ints(cls, values):
        for k, v in values.items():
            if k != "readout_os" and not isinstance(v, int):
                raise TypeError(f"Extra setting '{k}' must be an integer.")
        return values
    
    def get(self, key: str) -> int:
        if hasattr(self, key):
            return getattr(self, key)
        raise KeyError(f"Key '{key}' not found in settings.")

    def set(self, key: str, value: int):
        if not isinstance(value, int):
            raise TypeError(f"Value for '{key}' must be an integer.")
        setattr(self, key, value)


### infos section
class Info(BaseModel):
    description: str = "default_description"
    slices: int = 9999
    fov: int = 9999 # in [mm]
    dz: int = 5 # slice thickness in mm
    is3D: bool = False
    pelines: int = 9999
    seqstring: str = "default_seqstring"
    reconstruction: str = "default_reconstruction"


### instructions section
step_subclass_registry = {}


class Step(BaseModel):
    action: str
    
    def __init_subclass__(cls, **kwargs: Any) -> None:
        super().__init_subclass__(**kwargs)
        step_subclass_registry[cls.__name__] = cls    

    class Config:
        extra = "allow" 


class HasSteps():
    def __init__(self, **kwargs):
        for index in range(len(kwargs['steps'])):
            current_step = kwargs['steps'][index]
            if isinstance(current_step, dict):
                item_step_keys = sorted(current_step.keys())
                for name, cls in step_subclass_registry.items():
                    registery_step_keys = sorted(cls.model_fields.keys())
                    if item_step_keys == registery_step_keys:
                        try:
                            current_step = cls(**current_step)
                        except: 
                            pass
                        break
                else:
                    raise Exception(f"Unknown step action \"{current_step['action']}\"")
                kwargs['steps'][index] = current_step
        super().__init__(**kwargs)

    # def __init__(self, **kwargs):
    #     if(len(kwargs['steps'])!=0):
    #         for index in range(len(kwargs['steps'])):
    #             current_step = kwargs['steps'][index]
    #             if isinstance(current_step, dict):
    #                 item_step_keys = sorted(current_step.keys())
    #                 for name, cls in step_subclass_registry.items():
    #                     registery_step_keys = sorted(cls.model_fields.keys())
    #                     if item_step_keys == registery_step_keys:
    #                         try:
    #                             current_step = cls(**current_step)
    #                         except: 
    #                             pass
    #                         break
    #                 else:
    #                     raise Exception(f"Unknown step action \"{current_step['action']}\"")
    #                 kwargs['steps'][index] = current_step
    #         super().__init__(**kwargs)

class Instruction(HasSteps,BaseModel):
    print_counter: Optional[str] = "off"
    print_message: str = "default_print_message"

    steps: List[SerializeAsAny[Step]]


class Loop(HasSteps,Step):
    action: Literal["loop"] = "loop"
    counter: int = 9999
    range: int = 9999
    steps: List[SerializeAsAny[Step]] = []


class RunBlock(Step):
    action: Literal["run_block"] = "run_block"
    block: str = "default_block"


class Calc(Step):
    action: Literal["calc"] = "calc"
    type: str = "default_type"
    float: Decimal = 9.999
    increment: int = 9999


class Init(Step):
    action: Literal["init"] = "init"
    gradients: str = "default_gradients"


class Sync(Step):
    action: Literal["sync"] = "sync"
    object: str = "default_object"
    time: Union[int, EquationRef] = 9999


class Grad(Step):
    action: Literal["grad"] = "grad"
    axis: str = "default_axis"
    object: str = "default_object"
    time: Union[int, EquationRef] = 9999


class GradWithAmplitude(Grad):
    amplitude: Union[str, EquationRef] = "default_amplitude"


class EquationRef(BaseModel):
    type: str = "equation"
    equation: str = "default_equation"


class Rf(Step):
    action: Literal["rf"] = "rf"
    object: str = "default_object"
    time: Union[int, EquationRef] = 9999
    added_phase: AddedPhase


class Adc(Step):
    action: Literal["adc"] = "adc"
    object: str = "default_object"
    time: Union[int, EquationRef] = 9999
    frequency: int = 9999
    phase: int = 9999
    added_phase: AddedPhase  
    mdh: dict[str, MdhOption]


class AddedPhase(BaseModel):
    type: str = "default_added_phase_type"
    float: Decimal = 9.999


class MdhOption(BaseModel):
    type: str = "default_mdh_type"
    counter: int = 9999
    target: Optional[int] = None

class Mark(Step):
    action: Literal["mark"] = "mark"
    time: Union[int, EquationRef] = 9999


class Submit(Step):
    action: Literal["submit"] = "submit"


### objects section
object_subclass_registry = {}


class Object(BaseModel):
    type: str = "default_type"
    duration: int = 9999
    
    def __init_subclass__(cls, **kwargs: Any) -> None:
        super().__init_subclass__(**kwargs)
        object_subclass_registry[cls.__name__] = cls    

    class Config:
        extra = "allow" 


class HasObjects():
    def __init__(self, **kwargs):
        listed_objects = list(kwargs['objects'])
        for index in range(len(listed_objects)):
            current_object = listed_objects[index]
            if isinstance(current_object, dict):
                item_object_keys = sorted(current_object.keys())
                for name, cls in object_subclass_registry.items():
                    registery_step_keys = sorted(cls.model_fields.keys())
                    if item_object_keys == registery_step_keys:
                        try:
                            current_object = cls(**current_object)
                        except: 
                            pass
                        break
                else:
                    raise Exception(f"Unknown step action \"{current_object['objects']}\"")
                listed_objects[index] = current_object
        super().__init__(**kwargs)


class RfExcitation(Object):
    type: Literal["rf"] = "rf"
    array: str = "default_array"
    duration: int = 9999
    initial_phase: int = 9999
    thickness: int = 9999
    flipangle: int = 9999
    purpose: str = "default_purpose"


class GradientObject(Object):
    type: Literal["grad"] = "grad"
    array: str = "default_array"
    duration: int = 9999
    tail: int = 9999
    amplitude: float = 9.999


class AdcReadout(Object):
    type: Literal["adc"] = "adc"
    samples: int = 9999
    dwelltime: int = 9999


class Ttl(Object):
    type: Literal["sync"] = "sync"
    duration: int = 9999
    event: str = "default_event"


### arrays section
class GradientTemplate(BaseModel):
    encoding: str = "default_encoding"
    type: str = "default_type"
    size: int = 9999
    data: List[float] = []


### equations section
class Equation(BaseModel):
    equation: str = "default_equation"


### definition of entire SDL model
class PulseSequence(HasObjects, BaseModel):
    file: File
    settings: Settings
    infos: Info
    instructions: dict[str, Instruction]
    objects: dict[str, SerializeAsAny[Object]]
    arrays: dict[str, GradientTemplate]
    equations: dict[str, Equation]

### definition of entire SDL model
class PulseSequence(BaseModel):
    file: File = File()
    settings: Settings = Settings()
    infos: Info = Info()
    instructions: dict[str, Instruction] = {}
    objects: dict[str, SerializeAsAny[Object]] = {}
    arrays: dict[str, GradientTemplate] = {}
    equations: dict[str, Equation] = {}
  