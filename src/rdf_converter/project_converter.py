from pathlib import Path
import xml.etree.ElementTree as ET
from .utils.logging import get_logger

logger = get_logger(__name__)

class ProjectConverter:
    def __init__(self, xml_path: str):
        self.xml_path = Path(xml_path)

    def parse(self):
        logger.info(f"Parsing project XML: {self.xml_path}")
        self.tree = ET.parse(self.xml_path)
        self.root = self.tree.getroot()

    def to_turtle(self, out_ttl: str):
        Path(out_ttl).parent.mkdir(parents=True, exist_ok=True)
        with open(out_ttl, "w", encoding="utf-8") as w:
            w.write("@prefix jpost: <http://jpost.org/ontology/> .\n")
            w.write("@prefix : <http://jpost.org/resource/> .\n\n")
            # Minimal: just write project id/name if present
            proj_id = self.root.attrib.get("id", "PROJECT")
            w.write(f":{proj_id} a jpost:Project .\n")
