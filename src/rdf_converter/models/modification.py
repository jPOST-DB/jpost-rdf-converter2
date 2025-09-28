from __future__ import annotations

from dataclasses import dataclass, field
from ..utils.string_tool import is_not_empty


@dataclass
class Modification:
    title: str | None = None
    unimod: str | None = None
    site: str | None = None
    clazz: str | None = None

    def __init__(self):
        self.title = None
        self.unimod = None
        self.site = None
        self.clazz = None

    def get_title(self) -> str | None:
        return self.title
    
    def set_title(self, title: str) -> None:
        self.title = title

    def get_unimod(self) -> str | None:
        return self.unimod
    
    def set_unimod(self, unimod: str) -> None:
        self.unimod = unimod

    def get_site(self) -> str | None:
        return self.site
    
    def set_site(self, site: str) -> None:
        self.site = site

    def get_class(self) -> str | None:
        return self.clazz
    
    def set_class(self, clazz: str) -> None:
        self.clazz = clazz

    def to_ttl(self, f) -> None:
        class_id = 'JPO_034'

        if self.clazz == 'Post-translational':
            class_id = "JPO_021"
        elif self.clazz == 'Co-translational':
            class_id = "JPO_022"
        elif self.clazz == 'Pre-translational':
            class_id = "JPO_024"
        elif self.clazz.startswith("Chemical derivative"):
            class_id = "JPO_025"
        elif self.clazz == "Artefact":
            class_id = "JPO_026"
        elif self.clazz == "N-linked glycosylation":
            class_id = "JPO_027"
        elif self.clazz == "O-linked glycosylation":
            class_id = "JPO_028"
        elif self.clazz == "Other glycosylation":
            class_id = "JPO_029"
        elif self.clazz == "Synth. pep. protect. gp.":
            class_id = "JPO_030"
        elif self.clazz == "Isotopic label":
            class_id = "JPO_031"
        elif self.clazz == "Non-standard residue":
            class_id = "JPO_032"
        elif self.clazz == "Multiple":
            class_id = "JPO_033"
        elif self.clazz == "AA substitution":
            class_id = "JPO_035"
        elif self.clazz == "Cross-link":
            class_id = "JPO_036"
        elif self.clazz == "CID cleavable cross-link":
            class_id = "JPO_037"
        elif self.clazz == "Photo cleavable cross-link":
            class_id = "JPO_038"
        elif self.clazz == "Other cleavable cross-link":
            class_id = "JPO_039"

        f.write(f'        a unimod:UNIMOD_{self.unimod} ;\n')
        if is_not_empty(self.site):
            f.write(f'        jpost:modificationSite "{self.site}" ;\n')
        f.write(f'        jpost:modificationClass jpost:{class_id}\n')
