import json,tempfile,unittest
from pathlib import Path
from fetch_ibex import verified_cache
from experiment import identity,sha

class VerifiedCacheTests(unittest.TestCase):
    def test_corruption_or_incomplete_attempt_requires_refetch(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p=root/'method/donor';p.mkdir(parents=True)
            result=p/'calls.tsv';result.write_text('original calls\n');config={'sample':'S'}
            marker=p/'COMPLETE';marker.write_text(json.dumps({'configuration_sha256':identity(config),'outputs':{'calls.tsv':sha(result)}}))
            manifest=p/'manifest.json';manifest.write_text(json.dumps({'status':'complete','configuration':config}))
            self.assertEqual(verified_cache(root),{'method/donor':sha(marker)})
            result.write_text('changed calls\n');self.assertEqual(verified_cache(root),{})
            result.write_text('original calls\n');manifest.write_text(json.dumps({'status':'running','configuration':config}))
            self.assertEqual(verified_cache(root),{})

if __name__=='__main__':unittest.main()
