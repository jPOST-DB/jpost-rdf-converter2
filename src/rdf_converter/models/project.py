from __future__ import annotations

from dataclasses import dataclass, field
from .dataset import DataSet

from ..utils.string_tool import is_not_empty


@dataclass
class Contributor:
    name: str | None = None
    organization: str | None = None
    investigator: str | None = None

    def __init__(self):
        self.name = None
        self.organization = None
        self.investigator = None

    def get_name(self) -> str | None:
        return self.name
    
    def set_name(self, name: str) -> None:
        self.name = name

    def get_organization(self) -> str | None:
        return self.organization
    
    def set_organization(self, organization: str) -> None:
        self.organization = organization

    def get_investigator(self) -> str | None:
        return self.investigator
    
    def set_investigator(self, investigator: str) -> None:
        self.investigator = investigator




@dataclass
class Project:
    id: str
    dataset: DataSet | None = None
    pxd: str | None = None
    title: str | None = None
    description: str | None = None
    date_submitted: str | None = None
    date: str | None = None
    contributors: list[Contributor] = field(default_factory=list)

    def __init__(self, id: str):
        self.id = id
        self.dataset = None
        self.contributors = []

    def get_id(self) -> str:
        return self.id
    
    def set_id(self, id: str) -> None:
        self.id = id
    
    def get_project_number(self) -> int:
        string = self.id.replace('JPST', '')
        return int(string)
    
    def get_title(self) -> str | None:
        return self.title

    def set_title(self, title: str) -> None:
        self.title = title

    def get_pxd(self) -> str | None:
        return self.pxd
    
    def set_pxd(self, pxd: str) -> None:
        self.pxd = pxd

    def get_description(self) -> str | None:
        return self.description
    
    def set_description(self, description: str) -> None:
        self.description = description

    def get_dataset(self) -> DataSet | None:
        return self.dataset
    
    def set_dataset(self, dataset: DataSet) -> None:
        self.dataset = dataset

    def get_date_submitted(self) -> str | None:
        return self.date_submitted
    
    def set_date_submitted(self, date_submitted: str) -> None:
        self.date_submitted = date_submitted

    def get_date(self) -> str | None:
        return self.date
    
    def set_date(self, date: str) -> None:
        self.date = date

    def get_contributors(self) -> list[Contributor]:
        return self.contributors

    def to_ttl(self, f) -> None:
        f.write(f':{self.get_id()}\n')
        f.write(f'    dct:title "{self.get_title()}" ;\n')
        f.write(f'    dct:identifier "{self.get_id()}" ;\n')
        f.write(f'    rdfs:label "{self.get_id()}" ;\n')
        f.write(f'    rdfs:seeAlso pxd:{self.get_pxd()} ;\n')
        f.write(f'    rdfs:seeAlso jrepo:{self.get_id()} ;\n')
        f.write(f'    dct:description "{self.get_description()}" ;\n')
        f.write(f'    dct:dateSubmitted "{self.get_date_submitted()}"^^xsd:date ;\n')
        f.write(f'    dct:date "{self.get_date()}"^^xsd:date .\n\n')

        for contributor in self.get_contributors():
            name = contributor.get_name()
            organization = contributor.get_organization()
            investigator = contributor.get_investigator()

            f.write('    dct:contributor [\n')
            f.write('        a obo:MS_1002037 ;\n')
            f.write('        a foaf:Person ;\n')

            if is_not_empty(organization):
                f.write(f'        vcard:organization-name "{organization}" ;\n')

            if is_not_empty(name):
                f.write(f'        foaf:name "{name}" ;\n')

            f.write('    ] ;\n')

            if is_not_empty(investigator):
                f.write('    dct:contributor [\n')
                f.write('        a obo:MS_1002332 ;\n')
                f.write('        a foaf:Person ;\n')
                f.write(f'        foaf:name "{investigator}" ;\n')
                f.write('    ] ;\n')

            f.write('    a jpost:Project .\n\n')


    @staticmethod
    def read_project(xml_path: str) -> Project:
        import xml.etree.ElementTree as ET

        tree = ET.parse(xml_path)
        root = tree.getroot()

        node = root.find('Project')
        id = node.get('id')
        pxid = node.get('pxid')
        created_date = node.get('createdDate')

        title = node.find('Title').text
        description = node.find('Description').text
        date = node.find('AnnouncementDate').text

        project = Project(id)
        project.set_pxd(pxid)
        project.set_title(title)
        project.set_description(description)
        project.set_date_submitted(created_date)
        project.set_date(date)

        contributors = node.findall('Contributor')
        for contrib in contributors:
            contributor = Contributor()
            name = contrib.find('Name')
            if name is not None and name.text:
                contributor.set_name(name.text)
            organization = contrib.find('Affiliation')
            if organization is not None and organization.text:
                contributor.set_organization(organization.text)
            investigator = contrib.find('PrincipalInvestigator')
            if investigator is not None and investigator.text:
                contributor.set_investigator(investigator.text)

            if is_not_empty(contributor.get_name()):
                project.get_contributors().append(contributor)

        return project