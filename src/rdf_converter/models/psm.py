from __future__ import annotations

from dataclasses import dataclass, field
from tokenize import String
from typing import ClassVar

from typing import TYPE_CHECKING

from rdf_converter.utils.logging import get_logger
if TYPE_CHECKING:
    from .dataset import DataSet
    from .peptide import Peptide
    from .enzyme import Enzyme
    from .rawdata_list import RawDataList, RawData


from .spectrum import Spectrum
from .modification import Modification
from ..utils.string_tool import is_not_empty

logger = get_logger(__name__)


@dataclass
class PsmModification:
    modification : Modification | None = None
    site: str | None = None
    position: str | None = None

    def __init__(self, modification: Modification | None, site: str | None, position: str | None):
        self.modification = modification
        self.site = site
        self.position = position

    def get_modification(self) -> Modification | None:
        return self.modification
    
    def get_site(self) -> str | None:
        return self.site
    
    def get_position(self) -> str | None:
        return self.position


@dataclass
class Psm:
    dataset: DataSet
    peptide: Peptide | None = None
    sequence: str | None = None
    id: str | None = None
    title: str | None = None
    properties: dict = field(default_factory=dict)
    obs_mz: str | None = None
    calc_mz: str | None = None
    charge: str | None = None
    jpost_score: str | None = None
    mod: str | None = None
    mod_detail: str | None = None
    rt: str | None = None
    score_map: dict = field(default_factory=dict)
    fdr: float = 0.0
    raw_file: str | None = None
    scan: str | None = None
    phospho_confirmed: str | None = None
    phospho_ambiguous: list[str] = field(default_factory=list)
    representative: bool = False
    spectrum: Spectrum | None = None
    modifications: list[PsmModification] = field(default_factory=list)

    counter: ClassVar[int] = 0

    def __init__(self, dataset: DataSet):
        Psm.counter += 1

        self.dataset = dataset
        self.id = f'PSM{dataset.get_number()}_{Psm.counter}'
        self.properties = {}
        self.representative = False
        self.score_map = {}
        self.phospho_ambiguous = []
        self.modifications = []


    def get_dataset(self) -> DataSet:
        return self.dataset
    
    def get_peptide(self) -> Peptide | None:
        return self.peptide
    
    def set_peptide(self, peptide: Peptide) -> None:
        self.peptide = peptide

    def get_sequence(self) -> str | None:
        return self.sequence
    
    def set_sequence(self, sequence: str) -> None:
        self.sequence = sequence

    def get_id(self) -> str | None:
        return self.id
    
    def get_title(self) -> str | None:
        return self.title
    
    def set_title(self, title: str) -> None:
        self.title = title
        self.set_properties(title)

    def get_obs_mz(self) -> str | None:
        return self.obs_mz
    
    def set_obs_mz(self, obs_mz: str) -> None:
        self.obs_mz = obs_mz

    def get_calc_mz(self) -> str | None:
        return self.calc_mz
    
    def set_calc_mz(self, calc_mz: str) -> None:
        self.calc_mz = calc_mz

    def get_charge(self) -> str | None:
        return self.charge

    def set_charge(self, charge: str) -> None:
        self.charge = charge

    def get_jpost_score(self) -> str | None:
        return self.jpost_score
    
    def set_jpost_score(self, jpost_score: str) -> None:
        self.jpost_score = jpost_score

    def get_mod(self) -> str | None:
        return self.mod
    
    def set_mod(self, mod: str | None) -> None:
        self.mod = mod

    def get_mod_detail(self) -> str | None:
        return self.mod_detail
    
    def set_mod_detail(self, mod_detail: str | None) -> None:
        self.mod_detail = mod_detail

    def get_rt(self) -> str | None:
        return self.rt

    def set_rt(self, rt: str | None) -> None:
        self.rt = rt

    def get_score_map(self) -> dict:
        return self.score_map

    def set_score_map(self, score_map: dict) -> None:
        self.score_map = score_map

    def get_property(self, key: str) -> str | None:
        return self.properties.get(key)

    def set_property(self, key: str, value: str) -> None:
        self.properties[key] = value

    def get_fdr(self) -> float:
        return self.fdr
    
    def set_fdr(self, fdr: float) -> None:
        self.fdr = fdr

    def get_raw_file(self) -> str | None:
        return self.raw_file
    
    def set_raw_file(self, raw_file: str) -> None:
        self.raw_file = raw_file

    def get_scan(self) -> str | None:
        return self.scan
    
    def set_scan(self, scan: str) -> None:
        self.scan = scan

    def is_representative(self) -> bool:
        return self.representative
    
    def set_representative(self, representative: bool) -> None:
        self.representative = representative

    def get_phospho_confirmed(self) -> str | None:
        return self.phospho_confirmed

    def set_phospho_confirmed(self, phospho_confirmed: str) -> None:
        self.phospho_confirmed = phospho_confirmed

    def get_phospho_ambiguous(self) -> list[str]:
        return self.phospho_ambiguous
    
    def add_phospho_ambiguous(self, phospho_ambiguous: str) -> None:
        self.phospho_ambiguous.append(phospho_ambiguous)

    def get_spectrum(self) -> Spectrum | None:
        return self.spectrum

    def set_spectrum(self, spectrum: Spectrum) -> None:
        self.spectrum = spectrum

    def get_modifications(self) -> list[PsmModification]:
        return self.modifications

    def set_properties(self, title: str) -> None:
        items = title.split(',')
        for item in items:
            if ':' in item:
                index = item.index(':')
                if index > 0:
                    key = item[:index].strip()
                    value = item[index + 1:].strip()
                    self.set_property(key, value)


    def to_ttl(self, f) -> list[str]:
        self.modifications = []
        f.write(f':{self.id} \n')
        f.write(f'    dct:identifier "{self.id}" ;\n')

        if self.is_representative():
            f.write(f'    jpost:representativePsm 1 ;\n')
        
        if self.spectrum is not None:
            f.write(f'    jpost:hasSpectrum bid:{self.spectrum.get_id()} ;\n')

        f.write('    sio:SIO_000216 [\n')
        f.write('        a jpost:UniScore ;\n')
        f.write(f'        sio:SIO_000300 {self.get_jpost_score()} ;\n')
        f.write('    ] ;\n')

        modifications = []
        mod_set = set()
        mod_set21 = set()
        not_found = []

        dataset = self.get_dataset()
        enzyme = dataset.get_enzyme()
        for mod in enzyme.get_fixed_mods():
            modifications.append(mod)
        for mod in enzyme.get_variable_mods():
            modifications.append(mod)

        mod_map = {}
        tokens = self.get_mod().split(';')
        for token in tokens:
            mod = token
            spece_index = token.find(' ')
            if spece_index > 0:
                try:
                    int(token[:spece_index])
                    mod = token[spece_index + 1:]
                except ValueError:
                    pass

            start_index = token.find('(')
            end_index = token.find(')')
            mod_sites = token[start_index + 1: end_index]

            if mod_sites.find('N-term') >= 0:
                mod_map['N-term'] = mod
                mod_sites = mod_sites.replace('N-term', '')
            if mod_sites.find('C-term') >= 0:   
                mod_map['C-term'] = mod
                mod_sites = mod_sites.replace('C-term', '')
            while len(mod_sites) > 0:
                mod_site = mod_sites[0]
                mod_map[mod_site] = mod
                mod_sites = mod_sites[1:]

        tokens = self.mod_detail.split(',')
        for token in tokens:
            details = token
            index = details.find(':')
            position = None
            site = None
            mod = None
            
            if index >= 0:
                position = details[index + 1]
                index2 = details.find('@')
                if index2 >= 0 and index2 < index:
                    site = details[index2 + 1: index]
                    mod = mod_map.get(site)

            mod_list = []
            if mod is not None:
                if mod.find('(STY)') >= 0:
                    mod_list.append(mod.replace('(STY)', '(T)'))
                else:
                    mod_list.append(mod)

            for mod_element in mod_list:
                mod_info = f'Mod:{mod_element}:{site}:{position}'

                modification = None
                for tmp in modifications:
                    if tmp.get_title() == mod_element:
                        modification = tmp
                        if site is None:
                            site = tmp.get_site()

                if mod_info not in mod_set:
                    mod_set.add(mod_info)

                    f.write('    jpost:hasModification [\n')
                    if modification is not None:
                        if modification.get_unimod() == '21':
                            mod_info21 = f'Site:{site}, Position: {position}'
                            mod_set21.add(mod_info21)
                        f.write(f'        a unimod:UNIMOD_{modification.get_unimod()}\n')
                        psm_modification = PsmModification(modification, site, position)
                        self.modifications.append(psm_modification)
                    else:
                        if mod_element not in not_found:
                            not_found.append(mod_element)                            
                        f.write(f'        rdfs:label "{mod_element} " ;\n')

                    if site is not None:
                        f.write(f'        jpost:modificationSite "{site}" ;\n')
                    
                    if position is not None:
                        f.write('        faldo:location [\n')
                        f.write('            a faldo:ExactPosition ;\n')
                        f.write(f'            faldo:reference :{self.peptide.get_id()} ;\n')
                        f.write(f'            faldo:position {position} ;\n')
                        f.write('        ] \n')
                    f.write('    ] ;\n')

        f.write('    sio:SIO_000216 [\n')
        f.write('        a jpost:ExperimentalMassToCharge ;\n')
        f.write('        sio:SIO_000221 obo:MS_1000040 ;\n')
        f.write(f'        sio:SIO_000300 {self.get_obs_mz()};\n')
        f.write('   ] ;\n')

        f.write('    sio:SIO_000216 [\n')
        f.write('        a jpost:CalculatedMassToCharge ;\n')
        f.write('        sio:SIO_000221 obo:MS_1000040 ;\n')
        f.write(f'        sio:SIO_000300 {self.get_calc_mz()};\n')
        f.write('   ] ;\n')

        f.write('    sio:SIO_000216 [\n')
        f.write('        sio:SIO_000221 obo:MS_1000041 ;\n')
        f.write(f'        sio:SIO_000300 {self.get_charge()};\n')
        f.write('   ] ;\n')

        f.write('    sio:SIO_000216 [\n')
        f.write('        sio:SIO_000221 obo:MS_1000894 ;\n')
        f.write(f'        sio:SIO_000300 {self.get_rt()};\n')
        f.write('   ] ;\n')

        score_map = self.get_score_map()
        ev = score_map.get('ev')
        score = None
        score_id = ''
        ev_id = ''

        if 'Mascot' in score_map:
            score = score_map['Mascot']
            score_id = 'MS_1001171'
            ev_id = 'MS_1001172'
        elif 'X!Tandem' in score_map:
            score = score_map['X!Tandem']
            score_id = 'MS_1001331'
            ev_id = 'MS_1001330'
        elif 'Coment' in score_map:
            score = score_map['Coment']
            score_id = 'MS_1002252'
            ev_id = 'MS_1002257'
        elif 'MaxQuant' in score_map:
            score = score_map['MaxQuant']
            score_id = 'MS_1002338'
            ev_id = 'MS_1001901'

        if score_id is not None:
            f.write('    sio:SIO_000216 [\n')
            f.write(f'        a obo:{score_id} ;\n')
            f.write(f'        sio:SIO_000300 {score}\n')
            f.write('   ] ;\n')

        if ev_id is not None:
            if ev is not None and ev_id != '':            
                f.write('    sio:SIO_000216 [\n')
                f.write(f'        a obo:{ev_id} ;\n')
                f.write(f'        sio:SIO_000300 {ev}\n')
                f.write('   ] ;\n')

        self.write_phospho(f, mod_set21)
        f.write('    a jpost:Psm .\n\n')
        
        return not_found


    def write_phospho(self, f, mod_set21) -> None:
        not_found = []
        confirmed = self.get_phospho_confirmed()
        if is_not_empty(confirmed):
            array = confirmed.split('/')
            for element in array:
                values = element.split(':')
                site = ''
                pos = ''
                if len(values) >= 2:
                    site = values[0]
                    pos = values[1]

                mod_info21 = f'Site:{site}, Position:{pos}'
                if is_not_empty(site) and is_not_empty(pos) and mod_info21 not in mod_set21:
                    f.write('    jpost:hasModification [\n')
                    f.write('        a jpost:Modification ;\n')
                    f.write('        a unimod:UNIMOD_21 ;\n')
                    f.write(f'        jpost:modificationSite "{site}" ;\n')
                    f.write('        faldo:location [\n')
                    f.write('            a faldo:ExactPosition ;\n')
                    f.write(f'            faldo:reference :{self.get_peptide().get_id()} ;\n')
                    f.write(f'            faldo:position {pos} ;\n')
                    f.write('        ]\n')
                    f.write('    ] ;\n')

        ambiguousList = self.get_phospho_ambiguous()
        for ambiguous in ambiguousList:
            if ambiguous is not None and ambiguous != '':
                left_right = ambiguous.split('|')
                left = ''
                right = ''
                if len(left_right) >= 2:
                    left = left_right[0]
                    right = left_right[1]
                else:
                    right = ''

                if left:
                    f.write('    jpost:hasModification [\n')
                    f.write('        a jpost:AmbiguousModification ;\n')
                    f.write('        a unimod:UNIMOD_21 ;\n')
                    if left.startswith("!"):
                        f.write('        jpost:hasCorrespondingConfirmedSite true ;\n')
                        left = left[1:].strip()
                    else:
                        f.write('        jpost:hasCorrespondingConfirmedSite false ;\n')
                    site = ''
                    pos = ''
                    site_pos = left.split(':')
                    if len(site_pos) >= 2:
                        site = site_pos[0]
                        pos = site_pos[1]
                        mod_info21 = f"Site:{site}, Position:{pos}"
                        mod_set21.add(mod_info21)
                        if is_not_empty(site) and is_not_empty(pos):
                            f.write('        jpost:detectedSiteBySearchEngine [\n')
                            f.write(f'            jpost:modificationSite "{site}" ;\n')
                            f.write('            faldo:location [\n')
                            f.write('                a faldo:ExactPosition ;\n')
                            f.write(f'                faldo:reference :{self.get_peptide().get_id()} ;\n')
                            f.write(f'                faldo:position {pos} ;\n')
                            f.write('            ]\n')
                            f.write('        ] ;\n')
                    if is_not_empty(right):
                        array = right.split('+')
                        f.write("        faldo:location [\n")
                        f.write("            a faldo:OneOfPosition ;\n")
                        array_counter = 0
                        for array_element in array:
                            array_counter += 1
                            site_pos = array_element.split(":")
                            site = ''

                    if is_not_empty(right):
                        array = right.split("+")
                        f.write("        faldo:location [\n")
                        f.write("            a faldo:OneOfPosition ;\n")
                        array_counter = 0
                        for array_element in array:
                            array_counter += 1
                            site_pos = array_element.split(":")
                            site = ''
                            pos = ''
                            if len(site_pos) >= 2:
                                site = site_pos[0]
                                pos = site_pos[1]
                            if is_not_empty(pos):
                                f.write("            faldo:possiblePosition [\n")
                                f.write("                a faldo:ExactPosition ;\n")
                                f.write(f"                faldo:reference :{self.get_peptide().get_id()} ;\n")
                                f.write(f"                faldo:position {pos} ;\n")
                                f.write(f"                 jpost:modificationSite \"{site}\" ;\n")
                                f.write("            ]\n")

                                if array_counter < len(array):
                                    f.write(" ;\n")
                                f.write("\n")
                    f.write("        ]\n")
                    f.write("    ] ;\n")




    @staticmethod
    def get_psms(peptides: list[Peptide]) -> tuple[list[Psm], list[Spectrum]]:
        psms = []
        spectra = []
        spectra_map = {}
        rawdata_list = peptides[0].get_dataset().get_rawdata_list()

        for peptide in peptides:
            for psm in peptide.get_psms():
                psms.append(psm)
                file = psm.get_raw_file()
                rawdata = None
                for tmp in rawdata_list.get_list():
                    if rawdata is None or tmp.get_file() == file:
                        rawdata = tmp
            
                spectrum_id = f'{rawdata.get_id().replace("RAW", "SPC")}_{psm.get_scan()}'
                spectrum = None
                if spectrum_id in spectra_map:
                    spectrum = spectra_map[spectrum_id]
                else:
                    spectrum = Spectrum(rawdata, psm.get_scan())
                    spectra_map[spectrum_id] = spectrum
                    spectra.append(spectrum)
                psm.set_spectrum(spectrum)

        return psms, spectra
    
    @staticmethod
    def save_modifications(f, psms: list[Psm]) -> None:
        headers = ['Peptide ID', 'Modification', 'Site', 'Position']
        f.write('\t'.join(headers) + '\n')

        lines = []
        for psm in psms:
            peptide = psm.get_peptide()
            for psm_modification in psm.get_modifications():
                modification = psm_modification.get_modification()
                site = psm_modification.get_site()
                position = psm_modification.get_position()
                if modification is not None:
                    row = f'{peptide.get_id()}\t{modification.get_title()}\t{site}\t{position}'
                    lines.append(row)

        for line in sorted(lines):
            f.write(f'{line}\n')
