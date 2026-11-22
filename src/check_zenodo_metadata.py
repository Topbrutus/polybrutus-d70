from pathlib import Path
import json
R=Path(__file__).resolve().parents[1]
c=json.loads((R/'results/certificate_d70.json').read_text('utf8'))
m=json.loads((R/'publication/zenodo_metadata_blueprint.json').read_text('utf8'))
assert c['status']=='PASS' and c['dimension']==m['scientific_dimension']==70
assert c['vertices']==m['norm_one_vertices']==75
assert c['pairs']==m['exact_labelled_squared_pairs']==2775
assert c['affine_rank']==c['centered_gram_rank']==70
assert c['eigenvalue_75_over_2_multiplicity']==m['model_eigenvalue_75_over_2_multiplicity']==64
assert c['eigenvalue_37_multiplicity']==m['model_eigenvalue_37_multiplicity']==4
assert c['parent_tag_commit']=='3c8d8de8a19d458927d01bdfdbc811ab4b6eee6f'
assert c['parent_doi']==m['scientific_parent_doi']=='10.5281/zenodo.24433306'
assert m['related_identifiers']==[{'identifier':c['parent_doi'],'relation':'isDerivedFrom','resource_type':'publication-preprint','scheme':'doi'}]
assert len(m['keywords'])==9 and len(set(m['keywords']))==9
assert 'affine Euclidean dimension seventy' in m['keywords'] and 'sixty-fourfold ideal eigenvalue' in m['keywords']
assert not set(m['keywords']).intersection({'affine Euclidean dimension sixty-nine','sixty-threefold ideal eigenvalue','affine Euclidean dimension seventy-one','sixty-fivefold ideal eigenvalue'})
assert len(m['references'])==2 and c['parent_doi'] in m['references'][0] and '10.5281/zenodo.24407952' in m['references'][1]
assert all(x in (R/'publication/POLYBRUTUS_70D_NOTE.md').read_text('utf8') for x in ('## References',c['parent_doi'],'10.5281/zenodo.24407952'))
print('D70 METADATA 9 KEYWORDS / 2 REFERENCES / DIRECT D69 DOI PASS')
