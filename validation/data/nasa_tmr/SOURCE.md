# NASA Turbulence Modeling Resource (TMR): 2DZP, 2DBC, 2DN00

- **Fetched:** 2026-10-29
- **Site:** https://tmbwg.github.io/turbmodels/ (the old `turbmodels.larc.nasa.gov` now 301-redirects to a nasa.gov landing page that links here)
- **Repository:** https://github.com/TMBWG/turbmodels, pinned to commit `20a39f549dbff988a6accc421ac84dafc4487695` (committed 2026-09-02T18:58:55Z)
- **License:** CC0-1.0 (repository license per the GitHub API). NASA work, public domain in the US.
- **Grids:** taken from the official grid zips on nasa.gov, reading only the needed zip members with HTTP Range requests:
  - https://www.nasa.gov/wp-content/uploads/2026/02/flatplategrids-grids.zip (16.4 MB)
  - https://www.nasa.gov/wp-content/uploads/2026/02/bumpgrids-grids.zip (166 MB)
  - https://www.nasa.gov/wp-content/uploads/2026/02/naca0012-grids.zip (110 MB)
  - The coarse levels in the zips are byte-identical to the same files in the GitHub repo. Checked by sha256 for 35x25, 89x41 and 113-33.
- **Reference results and inputs:** `https://raw.githubusercontent.com/TMBWG/turbmodels/<commit>/<path>`. Paths are in the table below: `FlatPlate/SST/...`, `FlatPlate_validation/...`, `Bump/SST|SSTV|SST2003/...` and `NACA0012_validation/...`.
- **Model equations:** `model_docs/sst.html` (TMR SST page, saved as HTML) and `model_docs/noteonrunning.html` (freestream-turbulence guidance).

## Layout

| dir | content |
|---|---|
| `2DZP_flatplate/grids` | Plot3D 2D grids (formatted, double precision): 35x25, 69x49, 137x97, 273x193 (gz), plus neutral map files (`.nmf`) |
| `2DZP_flatplate/sst` | CFL3D/FUN3D SST-Vm results: Cf(x=0.97) and CD grid convergence, Cf(x), u+/y+, k/omega/mut profiles, F1/F2, typical inputs; `cf_incomp_results_sstv.dat` is from **incompressible** codes (FUN3D-incomp, SC/Tetra, GHOST) |
| `2DBC_bump/grids` | 89x41, 177x81, 353x161 (gz) + nmf |
| `2DBC_bump/sst`, `sstv`, `sst2003` | CFL3D/FUN3D results for SSTm, SST-Vm and SST-2003m: force convergence (CL, CD, CDp, CDv), Cf convergence at x=0.6321975/0.75/0.8678025, Cf(x), Cp(x), profiles; `sst/cf_convergence_sst_{overflow,tau}.dat` add OVERFLOW and TAU |
| `2DN00_naca0012/grids` | 113-33, 225-65, 449-129 (gz), **897-257 (gz)**, nmf, surface definition `n0012points_superbig_clust_fix.dat` |
| `2DN00_naca0012/sst` | CFL3D SST and SST-V CL/CD/Cp/Cf at alpha = 0/10/15 on 897x257; typical CFL3D/FUN3D inputs |
| `2DN00_naca0012/experiment` | Ladson (tripped, Re=6e6, M=0.15), Gregory & O'Reilly, Abbott & von Doenhoff, McCroskey fit |

The 897x257 NACA grid goes one level past the "smallest 2-3" brief. I added it because every TMR SST reference number for 2DN00 was computed on it, and it is only 3.8 MB.

## Files
| file | bytes | sha256 |
|---|---:|---|
| `2DBC_bump/grids/bump_2levelsdown_353x161.nmf` | 1010 | `13685315a88529465a0f52d85f97cdd5a8e671972d9791130260e6b11c9fea07` |
| `2DBC_bump/grids/bump_2levelsdown_353x161.p2dfmt.gz` | 421624 | `e0b49f79da69ff2cdab2f0d7f8628af2a7823b6afb0e1d48d587af0186ba6689` |
| `2DBC_bump/grids/bump_3levelsdown_177x81.nmf` | 1010 | `cb58faf80a271522d31cc2d91153d5dd39db86b30975119f586c765ecb278050` |
| `2DBC_bump/grids/bump_3levelsdown_177x81.p2dfmt` | 697772 | `0f11cf2818516326be73b4523d84b7132f395df354b300aa1caf3971bff16902` |
| `2DBC_bump/grids/bump_4levelsdown_89x41.nmf` | 1010 | `91b7da7261d0b22891a362b8829124f53c673a754e7784612d3ebf5456d922d1` |
| `2DBC_bump/grids/bump_4levelsdown_89x41.p2dfmt` | 177623 | `0ce600d67ef2dced8c14ee936f3f68578767541deeafa718ce32cabefd7486ba` |
| `2DBC_bump/sst/bump_u_sst.dat` | 62880 | `78532273b8bb1d0ec25726d609393fd1145b47b5217a18d8bd5fa459b6ee48e0` |
| `2DBC_bump/sst/cf_bump_sst.dat` | 42357 | `a50ef4603a7cabf4f79c4b492fc86316dc310a27162560348629ad18a62a4bb7` |
| `2DBC_bump/sst/cf_convergence_sst.dat` | 1863 | `792ba3d94315bc35524ee8069417c3fa935d7a63f5a3d30b5d7accf634648989` |
| `2DBC_bump/sst/cf_convergence_sst_overflow.dat` | 1135 | `0a4eab5a64d2f084b1bc8a407ac901d21aa589753d1eec35c9ca13128169683b` |
| `2DBC_bump/sst/cf_convergence_sst_tau.dat` | 1021 | `05c0fcf30dd50af74979b72ab1669a863db8b0ca796ec1a8dfde745dc2fa2028` |
| `2DBC_bump/sst/cfl3d_typical.inp` | 4274 | `6a6ecb6fbda08a8531bac890ef264092e8ac64517c7c54c0fdfd62d53ce55d07` |
| `2DBC_bump/sst/cp_bump_sst.dat` | 42555 | `a449db9c7b7c54ab54888a6413494367a6f0dfe3c1f3b327a819b2ef51cc1f2e` |
| `2DBC_bump/sst/find_peak_mut_sst.dat` | 43679 | `2f6e5099b5044ca7fcbf712541de009c2060123e7e5647a87a8c9ade35b1468e` |
| `2DBC_bump/sst/force_convergence_sst.dat` | 1195 | `2ba8fd025b1f9998e407fbfc2514672b139cc760d9a1b2aa907a92c13786a264` |
| `2DBC_bump/sst/fun3d.info` | 240 | `1dc4811607fc6ad537fbe44d5d5dc9e5f14390ed9f0c8c12a0f72c225f69154f` |
| `2DBC_bump/sst/fun3d.nml_typical` | 1288 | `276c2785d8ae9ef3e1fc8cba62e8bfbd4fb99719c4e4a2a0cc4a3e5872b97b0c` |
| `2DBC_bump/sst/mut_0.75_sst.dat` | 60969 | `36fd4de96d7ad3ce7560c3fcf3843f6f02855d2741f63e830afbaaa80653414a` |
| `2DBC_bump/sst/sst_omega_k.dat` | 83341 | `408f91ef1b650a1d0b26fecb019199eb23383b3cb0c70b40325bd97175026eb3` |
| `2DBC_bump/sst2003/cf_bump_sst2003.dat` | 52612 | `a21194119a1adf9464bac5dfd027a088b7727b57cecf99b70fb6152c994c7a0d` |
| `2DBC_bump/sst2003/cf_convergence_sst2003.dat` | 1917 | `62dd94da029a9acf09884b24dbad211f09f00b89fcd31b086cca2dfb857e687c` |
| `2DBC_bump/sst2003/cfl3d_typical.inp` | 4275 | `88ed4034c179bff333a188c96e5f731e685eed21e5582f614c7730e1b33051d4` |
| `2DBC_bump/sst2003/cp_bump_sst2003.dat` | 90708 | `3453b6b1498a4cae56e8a56f9d05c879ec5a98a5865861df96bf2d68d9ef9c06` |
| `2DBC_bump/sst2003/force_convergence_sst2003.dat` | 1215 | `e8a860067263d6b82122c347859c0de8c9c8a7abd3e1286d7388b1e592860b92` |
| `2DBC_bump/sst2003/mut_0.75_sst2003.dat` | 70465 | `20aa6a7e5eb3f4c847bf5920525f7bf100f404f6e79cddd72f10e98ec83ae78c` |
| `2DBC_bump/sstv/bump_u_sstv.dat` | 62880 | `af868ef77a92898e3260f7da6d3dd72486d909dacc29800815551d8776a0acf1` |
| `2DBC_bump/sstv/cf_bump_sstv.dat` | 42357 | `77dc70c4861827007ac09819120bc0c1084b7660ff64bc52947abf5aed60d414` |
| `2DBC_bump/sstv/cf_convergence_sstv.dat` | 1928 | `5d32fbf03043eba666558628849c6133582c80199dfefdbd23b9308fa6c818a0` |
| `2DBC_bump/sstv/cfl3d_typical.inp` | 4264 | `0c675d4e4cd9bfedc2907afad90632326751d44d7f89cbbdfc11e3d74e2900c8` |
| `2DBC_bump/sstv/cp_bump_sstv.dat` | 42553 | `343921427e50d8fcc4be77f4e36bc8832d6009c34cc3d4c0bb8d81de6403544c` |
| `2DBC_bump/sstv/force_convergence_sstv.dat` | 1306 | `a94857862d753fae88a07a4c2e6658d7a80e23653f796857a59263b1830ef26f` |
| `2DBC_bump/sstv/mut_0.75_sstv.dat` | 60905 | `40eff8cb43a3725e6570695b55d97fe37a6f5263c26e44c6494f15dd8d4ca493` |
| `2DN00_naca0012/experiment/0012.abbottdata.cd.dat` | 569 | `21321fe22883672c344f3cd1369c69cdaf045ed8496cd0a976b741645c79e3d9` |
| `2DN00_naca0012/experiment/0012.abbottdata.cl.dat` | 696 | `3ee6a0c198f6402c9fc5fd8affc1115c1e87ff943b8b02775df86a7cfa3c55f3` |
| `2DN00_naca0012/experiment/0012.mccroskeydata.cl.dat` | 250 | `bf5848c235cedffef4e307f5cb52be1e7146bd5e480c202824a59d6d9754924e` |
| `2DN00_naca0012/experiment/CLCD_Ladson_expdata.dat` | 1343 | `78cd2f6aa4968e80f44cbf6c96f699bd9c6e45681d958ec528a10cf72ed23357` |
| `2DN00_naca0012/experiment/CL_Gregory_expdata.dat` | 308 | `0b7381a45412de50fcf0e0dddf67f1e8e35ba1a278bd4734fb7f16ad08df1f74` |
| `2DN00_naca0012/experiment/CP_Gregory_expdata.dat` | 1608 | `38d8474bc67661c3cef5e4cd9eeddb95d4a436d06abc1f58a8518066f510647c` |
| `2DN00_naca0012/experiment/CP_Ladson.dat` | 6331 | `547e022634a360cd7697813543348cd447ee0038cdbc8f00a5296c7822ec7af8` |
| `2DN00_naca0012/grids/n0012_113-33.nmf` | 957 | `597865103baf72a0c78d0e42a7e2d09688ea78b6a20f8fbcac2baa936d1445a1` |
| `2DN00_naca0012/grids/n0012_113-33.p2dfmt` | 181516 | `669fff259b59ac4df557ad3bfef6be0794ac5eacc0f45ec1ce9715f75dfd1654` |
| `2DN00_naca0012/grids/n0012_225-65.nmf` | 957 | `c6c32668f19bbd1a89510dfeb3a43c924137ae90afa7772eeb83875f9d9b34b0` |
| `2DN00_naca0012/grids/n0012_225-65.p2dfmt` | 711788 | `66817d6bbc9353cb0ca9867f958da1bcce6485f7d68ec7038362a341bdb0fa8c` |
| `2DN00_naca0012/grids/n0012_449-129.nmf` | 957 | `648e7b4c368db4479db3cdcb3d0209590da991ad00e59865398330cbc4b45880` |
| `2DN00_naca0012/grids/n0012_449-129.p2dfmt.gz` | 949473 | `6d60dfa6ef9c7544ee156276fcc5487d7eacddadbdf6e81cc6f10e78fb03eff8` |
| `2DN00_naca0012/grids/n0012_897-257.nmf` | 957 | `ec897e29363eb72546730a83048f142f8219d33397407b9f8cf474f94fe10365` |
| `2DN00_naca0012/grids/n0012_897-257.p2dfmt.gz` | 3776075 | `c603ea77017ae8b97af581237b81a8251fd2d3af1c49c8c6f2cd2687aafcb28e` |
| `2DN00_naca0012/grids/n0012points_superbig_clust_fix.dat` | 824246 | `74e147db4fa052be7ab79465a756c55b39b1786bfd7141a86257e89586728ae1` |
| `2DN00_naca0012/sst/fun3d.nml_typical_sst_2` | 1078 | `c5c40ffe4834e1485c840c7a594d96738441ebef96d0eeeb376e1fc232193bb8` |
| `2DN00_naca0012/sst/n0012_cfl3d_typical_sst.inp` | 6060 | `e3d1b4e7c9600d03f90e5e6b5cd7a10d01ad041f704e3d60de4cacd57cb375a4` |
| `2DN00_naca0012/sst/n0012cf_cfl3d_sst.dat` | 37962 | `d8001cd09d0c573de63475f8e0e02a2a989216133af74b84876eba5c328da8be` |
| `2DN00_naca0012/sst/n0012cf_cfl3d_sstv.dat` | 37964 | `7f0d7d4da9b2a977e620c29eed733e1e06745d291ab54d1888f7137886e05728` |
| `2DN00_naca0012/sst/n0012clcd_cfl3d_sst.dat` | 417 | `373dc142018816628b7b53ad278327db731162988ba35a166cd9fdc131a7fd7a` |
| `2DN00_naca0012/sst/n0012clcd_cfl3d_sstv.dat` | 419 | `5165603932574c1f03fc164fdf96ae71c57643cfdeea12a56c2d7c870cb8ec5b` |
| `2DN00_naca0012/sst/n0012cp_cfl3d_sst.dat` | 48141 | `378f9dcd058ffbe900081f213ee5c7e24dfaec0d6a9e02f366f96e8d76bf179d` |
| `2DN00_naca0012/sst/n0012cp_cfl3d_sstv.dat` | 48143 | `46af599e6c920f143c18aa4bf588228646277f4182c2d051d3f72ee9c113f8dc` |
| `2DZP_flatplate/grids/flatplate_clust2_1leveldown_273x193.nmf` | 927 | `a96ea21baa762e2d5a4f5d0447ae2ebaf2fb5323f1a0f20c36c32e5c0740f498` |
| `2DZP_flatplate/grids/flatplate_clust2_1leveldown_273x193.p2dfmt.gz` | 42020 | `40be4c0a92b75892189263a5c3efb63fa112030912480427e1e0787c6cb97955` |
| `2DZP_flatplate/grids/flatplate_clust2_2levelsdown_137x97.nmf` | 927 | `4e2d2e8d5a0076754c06a310f6b7790096a59d5c02319452f81234865f13887d` |
| `2DZP_flatplate/grids/flatplate_clust2_2levelsdown_137x97.p2dfmt` | 646770 | `55ae965f2151e5644e3c7c1b6fecea7ce0a16cf4c34eae6997f27683a4e2ac44` |
| `2DZP_flatplate/grids/flatplate_clust2_3levelsdown_69x49.nmf` | 927 | `9117aaa4a80c97aa0be8a558db4c3582b7ec89d0b0214c61158ce51fad248a64` |
| `2DZP_flatplate/grids/flatplate_clust2_3levelsdown_69x49.p2dfmt` | 164580 | `18c1d1c82cf91bdad89b2c2963129fc992fcd8b98fda4b3e79b7397ee52ba223` |
| `2DZP_flatplate/grids/flatplate_clust2_4levelsdown_35x25.nmf` | 927 | `cae0d171eb6e62843f8f954205d7edca688b1ddc1737145acdc6bf85b0e6dab6` |
| `2DZP_flatplate/grids/flatplate_clust2_4levelsdown_35x25.p2dfmt` | 42622 | `41c409f49afca8affaafd4a0a49c6ba934eca9007da502fe0c20afbd091500b0` |
| `2DZP_flatplate/sst/cf_convergence_fun3d_sstv_tri.dat` | 405 | `fc9392568cb0d5612552d5414af21e1bda98b3a999f6154b510159a884bcc004` |
| `2DZP_flatplate/sst/cf_convergence_sstv.dat` | 666 | `e5cde959873804d9b37a83640535309e18d5f1bc269b7bb6a46c460c40cd2a30` |
| `2DZP_flatplate/sst/cf_incomp_results_sstv.dat` | 986 | `65bbdbae974af5297db4387d67f2f9ac8d3c46750815d684e108a52e850647b9` |
| `2DZP_flatplate/sst/cf_plate_fun3d_sstv_tri.dat` | 14797 | `8a95c87fb3877d5a239f7ac06fd94870a3310d401ced832cdbb7d3ab29419b44` |
| `2DZP_flatplate/sst/cf_plate_sstv.dat` | 32816 | `a56bd2d16818456af015b13839f0a991b4b85bfa6c18242e7c331a2c7a4d40e8` |
| `2DZP_flatplate/sst/cfl3d_typical.inp` | 4539 | `ddebda5abf2fb31faa8650ccb6f753eaf377fe35512e26ec4e6883d5183d8c41` |
| `2DZP_flatplate/sst/drag_convergence_fun3d_sstv_tri.dat` | 343 | `f915bef096c37a4f135b6c77fad5bce10ccb6ef54e3b1c26f0191b9b12a547fa` |
| `2DZP_flatplate/sst/drag_convergence_sstv.dat` | 604 | `85fdc40fdacb447f954cac6704d6722cbdcb5ef5717acdfba4f76b1a7cba1f05` |
| `2DZP_flatplate/sst/find_peak_mut_sstv.dat` | 16895 | `5b1207c6da015ac35e9a2a868880be9f032f29cabbb6480d7185b2aef53f341e` |
| `2DZP_flatplate/sst/flatplate_u+y+_sstv.dat` | 23881 | `4a5bf0bcd8a3888bd313d67fd4a78f5cbde826fa0b00c62ac68d342197a534b2` |
| `2DZP_flatplate/sst/flatplate_u_sstv.dat` | 23936 | `4e1e17a24ffb64708ce805e978cf83a3b4391585e0ded7d14dd9be6e059b1d5f` |
| `2DZP_flatplate/sst/fun3d.nml_typical_2013` | 1303 | `3d67e465596822811ecb18c3dc559ce5287f84e924a52b90b7736ba145733b0c` |
| `2DZP_flatplate/sst/mut_0.97_sstv.dat` | 37737 | `9378362512c4e805772ac5a0d752dc92dbfe1daf17711013b896516aa95163c2` |
| `2DZP_flatplate/sst/sstv_f1f2_function.dat` | 31059 | `aaec944650609a32e12a18578a9a3473a669922a9adb85c2983cd255eb396096` |
| `2DZP_flatplate/sst/sstv_omega_k.dat` | 56059 | `caf7254a3c21195bef90a109ea98b8a7187bec959410a56387cd6c21410085bc` |
| `2DZP_flatplate/sst/u+y+theory.dat` | 380438 | `e4439ae7f484cf63599e053dfa03e6746d4c4519b47d1345644db4010afc6e69` |
| `model_docs/noteonrunning.html` | 35301 | `30979bfc7284dd3855a6299190ccf1b5f095ad7aa9ca776d30bb7f8f2b3a4fc1` |
| `model_docs/sst.html` | 35218 | `ccc4bc8d999444104e4e41f2c610dc9de7e7cb213d3e8b7bdb1997c4b8d6ee9f` |
