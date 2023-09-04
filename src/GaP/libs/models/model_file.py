
from pydantic import BaseModel

class FileModel(BaseModel):
    """ folder/filename
    """
    folder: str
    filename: str