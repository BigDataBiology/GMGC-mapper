import pytest
import sys
sys.path.append(r'../src')
from alignment import identity_coverage


KNOWN_RESULTS=[
    {
    'query_nt': "ATGATGATGATGTGA"
  , 'query_aa': "MMMM"
  , 'target_nt': "ATGATGATGATGTGA"
  , 'target_aa': "MMMM"
  , 'result': "EXACT"
    },
    {
        'query_nt': "ATGCATACTTATCACACCTACCATACCTAC"
        , 'query_aa': "MHTYHTYHTY"
        , 'target_nt': "ATGCACACTTACCATACGTATCATACATAC"
        , 'target_aa': "MHTYHTYHTY"
        , 'result': "SIMILAR"
    },
    {
        'query_nt': "ATGCATACTTATCACACCTACCATACCTAC"
        , 'query_aa': "MHTYHTYHTY"
        , 'target_nt': "ATGCACTACCCCCATTATCACACGTAT"
        , 'target_aa': "MHYPHYHTY"
        , 'result': "MATCH"
    },
    {
        'query_nt': "ATGCATACTTATCACACCTACCATACCTAC"
        , 'query_aa': "MHTYHTYHTY"
        , 'target_nt': "ATGGAACCTGAGCCAGAACCC"
        , 'target_aa': "MEPEPEP"
        , 'result': "NO MATCH"
    }
]


@pytest.mark.parametrize("t", KNOWN_RESULTS)
def test_alignment(t):
        assert identity_coverage(t['query_nt'],t['query_aa'],t['target_nt'],t['target_aa']) == t['result']



