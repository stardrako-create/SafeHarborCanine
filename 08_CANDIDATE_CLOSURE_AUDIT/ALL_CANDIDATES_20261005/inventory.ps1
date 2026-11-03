param([string]$Project)
$ErrorActionPreference='Stop'
$o=Join-Path $Project '08_CANDIDATE_CLOSURE_AUDIT\ALL_CANDIDATES_20261005'
$f=Join-Path $Project '06_v1.21.8\candidates_integrated_review_v1218.tsv'
$r=@(Import-Csv -LiteralPath $f -Delimiter "`t")
if($r.Count -ne 461){throw 'Unexpected regional universe'}
$gates=@($r[0].PSObject.Properties.Name | Where-Object {$_ -like 'veto_*'})
$entries=@(foreach($x in $r){$failed=@($gates | Where-Object {$x.$_ -eq 'True'}); [ordered]@{chrom=$x.chrom;start0=[int]$x.start;end0=[int]$x.end;genes="$($x.left_gene)/$($x.right_gene)";hard_veto=$x.hard_veto;failed_gates=$failed;missing_evidence=$x.missing_evidence;integrated_review_status=$x.integrated_review_status}})
$counts=@(foreach($g in $gates){[ordered]@{gate=$g;failed=@($r | Where-Object {$_.$g -eq 'True'}).Count}})
$a=Join-Path $Project '07_FINAL_CANDIDATES_2026-09-19\all16_decisions_and_evidence.tsv'
$b=Join-Path $Project '07_FINAL_CANDIDATES_2026-09-19\remaining48_no_new_biological_veto.tsv'
$reviewed=@(Import-Csv $a -Delimiter "`t");$remaining=@(Import-Csv $b -Delimiter "`t")
if($reviewed.Count -ne 16 -or $remaining.Count -ne 48){throw 'Unexpected window universe'}
$all=@($reviewed)+@($remaining)
$keys=@($all | ForEach-Object {"$($_.chrom):$($_.start)-$($_.end)"} | Sort-Object -Unique)
if($keys.Count -ne 64){throw 'Nonunique window identities'}
[ordered]@{regional_count=$r.Count;window_count=$all.Count;unique_windows=$keys.Count;regional_status=@($r | Group-Object integrated_review_status | Select-Object Name,Count);gate_counts=$counts;regions=$entries;window_region_counts=@($all | Group-Object genes | Select-Object Name,Count);source_hashes=@($f,$a,$b | Get-FileHash -Algorithm SHA256 | Select-Object Path,Hash);limitations=@('Inventory of existing evidence, not new biological validation','Gate counts overlap','Regional vetoes remain in force; subwindows require their own complete evaluation')} | ConvertTo-Json -Depth 9 | Set-Content (Join-Path $o 'inventory.json') -Encoding utf8
