import os


class SaveConfig:
    """Singleton holding the user-configurable base directory for all data saves."""

    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance.base_dir = r"D:\pumpprobedata"
        return cls._instance

    def set_base_dir(self, path: str):
        self.base_dir = os.path.normpath(path)

    def date_dir(self, timestamp) -> str:
        """Return base_dir/YYYY/MM/DD for a given datetime object."""
        return os.path.join(
            self.base_dir,
            timestamp.strftime("%Y"),
            timestamp.strftime("%m"),
            timestamp.strftime("%d"),
        )
