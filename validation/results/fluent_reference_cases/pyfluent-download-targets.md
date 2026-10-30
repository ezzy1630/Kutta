# PyFluent download targets

Snapshot ee3529f845f290573d63758e9f78778f1b8b00d9; searched all 108 documentation and example source files. Code below is quoted MIT-licensed source text, never executed; copyright and permission notice retained in `LICENSE.pyfluent`.

## doc/source/user_guide/events.rst
45:   >>> case_file_name = examples.download_file(
46:   >>>     "mixing_elbow.cas.h5", "pyfluent/mixing_elbow"
47:   >>> )
48:   >>> def on_case_loaded(session, event_info):
49:   >>>     on_case_loaded.loaded = True
50:
51:   >>> on_case_loaded.loaded = False
80:   >>> case_file_name = examples.download_file(
81:   >>>     "mixing_elbow.cas.h5", "pyfluent/mixing_elbow"
82:   >>> )
83:   >>> def on_case_loaded(session, event_info):
84:   >>>     on_case_loaded.loaded = True
85:
86:   >>> on_case_loaded.loaded = False

## doc/source/user_guide/fields/field_data.rst
21:   >>> case_path = download_file(file_name="exhaust_system.cas.h5", directory="pyfluent/exhaust_system")
22:   >>> data_path = download_file(file_name="exhaust_system.dat.h5", directory="pyfluent/exhaust_system")
23:   >>> solver_session.settings.file.read_case_data(file_name=case_path)
24:
25:   >>> field_data = solver_session.fields.field_data  # This creates an instance of the FieldData class.
26:
27: Simple requests
28: ---------------
335:   >>> import_file_name = examples.download_file(
336:   >>>     "mixing_elbow.pmdb", "pyfluent/mixing_elbow"
337:   >>> )
338:
339:   >>> # Launch Fluent in Meshing mode
340:   >>> meshing_session = pyfluent.Meshing.from_install()
341:

## doc/source/user_guide/fields/reduction.rst
27:     >>> case_path = download_file(file_name="exhaust_system.cas.h5", directory="pyfluent/exhaust_system")
28:     >>> data_path = download_file(file_name="exhaust_system.dat.h5", directory="pyfluent/exhaust_system")
29:     >>> solver1.settings.file.read_case_data(file_name=case_path)
30:
31:     >>> solver2 = pyfluent.Solver.from_install()
32:     >>> case_path = download_file("elbow1.cas.h5", "pyfluent/file_session")
33:     >>> data_path = download_file("elbow1.dat.h5", "pyfluent/file_session")
34:     >>> solver2.settings.file.read_case_data(file_name=case_path)
35:
36:     >>> solver_session =solver1
37:
38:
39: Functional usage

## doc/source/user_guide/fields/solution_data.rst
21:   >>> import_file_name = examples.download_file("mixing_elbow.cas.h5", "pyfluent/mixing_elbow")
22:   >>> examples.download_file("mixing_elbow.dat.h5", "pyfluent/mixing_elbow")
23:   >>> solver_session = pyfluent.Solver.from_install()
24:   >>> solver_session.settings.file.read_case_data(file_name=import_file_name)
25:
26:
27: The ``solution_variable_info`` and ``solution_variable_data`` objects are attributes of the ``solver_session.fields`` object:
28:

## doc/source/user_guide/file_transfer.rst
29:       >>> case_file = examples.download_file("mixing_elbow.cas.h5", "pyfluent/mixing_elbow")
30:
31:       >>> session = pyfluent.Solver.from_install()
32:       >>> session.settings.file.read_case(file_name=case_file)
33:       >>> session.settings.file.write_case(file_name="write_mixing_elbow.cas.h5")
34:
35:
54:       >>> case_file = examples.download_file("mixing_elbow.cas.h5", "pyfluent/mixing_elbow", return_without_path=False)
55:
56:       >>> file_transfer_service = ContainerFileTransferStrategy()
57:       >>> container_dict = {"mount_source": file_transfer_service.mount_source}
58:       >>> session = pyfluent.Solver.from_container(file_transfer_service=file_transfer_service, container_dict=container_dict)
59:       >>> session.settings.file.read_case(file_name=case_file)
60:       >>> session.settings.file.write_case(file_name="write_mixing_elbow.cas.h5")
81:       >>> mesh_file = examples.download_file("mixing_elbow.msh.h5", "pyfluent/mixing_elbow")
82:
83:       >>> session = pyfluent.Meshing.from_install(file_transfer_service=StandaloneFileTransferStrategy())
84:       >>> session.upload(file_name=mesh_file, remote_file_name="elbow.msh.h5")
85:       >>> session.meshing.File.ReadMesh(FileName="elbow.msh.h5")
86:       >>> session.meshing.File.WriteMesh(FileName="write_elbow.msh.h5")
87:       >>> session.download(file_name="write_elbow.msh.h5", local_directory="<local_path>")

## doc/source/user_guide/legacy/meshing_workflows.rst
23:     >>> import_file_name = examples.download_file('mixing_elbow.pmdb', 'pyfluent/mixing_elbow')
24:     >>> meshing_session = pyfluent.launch_fluent(
25:     >>>     mode=pyfluent.FluentMode.MESHING, precision=pyfluent.Precision.DOUBLE, processor_count=2
26:     >>> )
27:     >>> workflow = meshing_session.workflow
28:     >>> workflow.InitializeWorkflow(WorkflowType='Watertight Geometry')
29:     >>> tasks = workflow.TaskObject
142:     >>> import_file_name = examples.download_file(
143:     >>>     "exhaust_system.fmd", "pyfluent/exhaust_system"
144:     >>> )
145:     >>> meshing_session = pyfluent.launch_fluent(
146:     >>>     precision=pyfluent.Precision.DOUBLE,
147:     >>>     processor_count=2,
148:     >>>     mode=pyfluent.FluentMode.MESHING
614:     >>> import_file_name = examples.download_file('NACA0012.fmd', 'pyfluent/airfoils')
615:     >>> meshing_session = pyfluent.launch_fluent(
616:     >>>     mode=pyfluent.FluentMode.MESHING,
617:     >>>     precision=pyfluent.Precision.DOUBLE,
618:     >>>     processor_count=2
619:     >>> )
620:     >>> workflow = meshing_session.workflow

## doc/source/user_guide/meshing/new_meshing_workflows.rst
34:     import_file_name = examples.download_file('mixing_elbow.pmdb', 'pyfluent/mixing_elbow')
35:     meshing_session = pyfluent.Meshing.from_install(
36:         precision=pyfluent.Precision.DOUBLE, processor_count=2
37:     )
38:     watertight = pyfluent.WatertightMeshing(session=meshing_session)
39:     import_geometry = watertight.import_geometry
40:     import_geometry.file_name = import_file_name
136:     import_file_name = examples.download_file(
137:         "exhaust_system.fmd", "pyfluent/exhaust_system"
138:     )
139:     meshing_session = pyfluent.Meshing.from_install(precision=pyfluent.Precision.DOUBLE, processor_count=2)
140:
141:     fault_tolerant = FaultTolerantMeshing(session=meshing_session)
142:     fault_tolerant.parts.input_file_changed(
391:     import_file_name = examples.download_file('NACA0012.fmd', 'pyfluent/airfoils')
392:     meshing_session = pyfluent.Meshing.from_install(
393:         precision=pyfluent.Precision.DOUBLE, processor_count=2
394:     )
395:     two_dim_mesh = TwoDimensionalMeshing(session=meshing_session)
396:
397:     load_cad = two_dim_mesh.load_cad_geometry
531:     import_file_name = examples.download_file('mixing_elbow.pmdb', 'pyfluent/mixing_elbow')
532:     meshing_session = pyfluent.Meshing.from_install(
533:         precision=pyfluent.Precision.DOUBLE, processor_count=2
534:     )
535:     custom_workflow = CreateMeshingWorkflow(session=meshing_session)
536:
537: Insert first task
576:     saved_workflow_path = examples.download_file(
577:         "sample_watertight_workflow.wft", "pyfluent/meshing_workflows"
578:     )
579:     meshing_session = pyfluent.Meshing.from_install(
580:         precision=pyfluent.Precision.DOUBLE, processor_count=2
581:     )
582:     loaded_workflow = LoadMeshingWorkflow(session=meshing_session, file_path=saved_workflow_path)

## doc/source/user_guide/monitors.rst
26:   >>> import_case = examples.download_file(
27:   >>>     file_name="exhaust_system.cas.h5", directory="pyfluent/exhaust_system"
28:   >>> )
29:   >>> import_data = examples.download_file(
30:   >>>     file_name="exhaust_system.dat.h5", directory="pyfluent/exhaust_system"
31:   >>> )
32:   >>> solver_session.settings.file.read_case_data(file_name=import_case)
33:   >>> # check the active report plot monitors using the settings relevant object
34:   >>> sorted(solver_session.settings.solution.monitor.report_plots())
35:   >>> [

## doc/source/user_guide/offline/case_file.rst
22:   >>> case_file_name = examples.download_file("Static_Mixer_Parameters.cas.h5", "pyfluent/static_mixer")
23:   >>> reader = CaseFile(case_file_name=case_file_name)
24:   >>> reader.precision()
25:   2
26:   >>> reader.num_dimensions()
27:   3
28:   >>> {p.name: p.value for p in reader.input_parameters()}
68:       >>> case_file_name = examples.download_file("elbow1.cas.h5", "pyfluent/file_session")
69:       >>> reader = CaseFile(case_file_name=case_file_name)
70:       >>> reader.get_mesh().get_surface_ids()
71:       [3, 4, 5, 6, 7, 9]
72:       >>> reader.get_mesh().get_surface_names()
73:       ['wall',
74:        'symmetry',

## doc/source/user_guide/offline/data_file.rst
23:   >>> data_file_name = examples.download_file("elbow1.dat.h5", "pyfluent/file_session")
24:   >>> reader = DataFile(data_file_name=data_file_name, case_file_handle=CaseFile(case_file_name))
25:   >>> reader.case_file
26:   'elbow1.cas.h5'
27:   >>> reader.variables()
28:   {'flow-time': 0.0,
29:    'solvertime': 43.13454608435059,

## doc/source/user_guide/offline/expression_builder.rst
25:   >>> case_path = download_file("Static_Mixer_main.cas.h5", "pyfluent/static_mixer")
26:   >>> solver_session.settings.file.read_case(file_name=case_path)
27:
28:   >>> # Initialize the solution before creating and evaluating expressions
29:   >>> solver_session.solution.initialization.hybrid_initialize()
30:
31:   >>> # Create the expression builder

## doc/source/user_guide/offline/file_session.rst
36:   >>> case_file_name = examples.download_file("elbow1.cas.h5", "pyfluent/file_session")
37:   >>> data_file_name = examples.download_file("elbow1.dat.h5", "pyfluent/file_session")
38:   >>> file_session = FileSession(case_file_name, data_file_name)
39:
40:   >>> file_session.fields.field_data.scalar_fields.range(VariableCatalog.TEMPERATURE)
41:   [0.0, 313.1515948109515]
42:   >>> file_session.fields.field_data.surfaces()
43:   {'wall': {'surface_id': [3],
135:   >>> case_file_name = examples.download_file("mixing_elbow_mul_ph.cas.h5", "pyfluent/file_session")
136:   >>> data_file_name = examples.download_file("mixing_elbow_mul_ph.dat.h5", "pyfluent/file_session")
137:   >>> file_session = FileSession()
138:   >>> file_session.read_case(case_file_name)
139:   >>> file_session.read_data(data_file_name)
140:
141:   >>> file_session.fields.field_data.scalar_fields.range("phase-2:SV_P")
142:   [0.0, 1.5435200335871788e-11]

## doc/source/user_guide/physical_variables.rst
34:     case_path = examples.download_file("mixing_elbow.cas.h5", "pyfluent/mixing_elbow")
35:     solver_session.settings.file.read(file_type="case", file_name=case_path)
36:
37:     solver_session.settings.solution.initialization.hybrid_initialize()
38:
39:     temperature = VariableCatalog.TEMPERATURE
40:     locations = ["hot-inlet"]

## doc/source/user_guide/session/launching_ansys_fluent.rst
507:   >>> case_file_name = examples.download_file("mixing_elbow.cas.h5", "pyfluent/mixing_elbow")
508:   >>> solver_session.settings.file.read(file_name=case_file_name, file_type="case")
509:   >>> solver_session.exit()
510:
511:
512: Connect to a Fluent container running inside WSL from a Windows host
513: --------------------------------------------------------------------

## doc/source/user_guide/session/session.rst
21:   >>> case_file_name = download_file("mixing_elbow.cas.h5", "pyfluent/mixing_elbow")
22:   >>> solver_session = pyfluent.Solver.from_install(case_file_name=case_file_name)
23:
24: .. note::
25:
26:   You can also use the lower-level :func:`launch_fluent() <ansys.fluent.core.execution.launcher.launcher.launch_fluent>`
27:   and :func:`connect_to_fluent() <ansys.fluent.core.execution.launcher.launcher.connect_to_fluent>` functions,
134:   >>> import_file_name = examples.download_file('mixing_elbow.pmdb', 'pyfluent/mixing_elbow')
135:   >>> import_geometry = watertight.import_geometry
136:   >>> import_geometry.file_name.set_state(import_file_name)
137:   >>> length_unit = import_geometry.length_unit
138:   >>> length_unit.get_state()
139:    "mm"
140:   >>> length_unit.allowed_values()
205:   >>> case_file = download_file("mixing_elbow.cas.h5", "pyfluent/mixing_elbow")
206:   >>> with using(solver_session_1):
207:   ...     ReadCase()(file_name=case_file)
208:   ...     Viscous().model.set_state("laminar")
209:   ...     print(Viscous().model())
210:   laminar
211:   >>> with using(solver_session_2):

## doc/source/user_guide/solver_settings/set_up/boundary_conditions.rst
13:     >>> file_name = examples.download_file("mixing_elbow.cas.h5", "pyfluent/mixing_elbow")
14:     >>> solver_session = pyfluent.Solver.from_install()
15:     >>> solver_session.settings.file.read_case(file_name=file_name)
16:     >>> cold_inlet = pyfluent.solver.VelocityInlet(settings_source=solver_session, name="cold-inlet")
17:     >>> cold_inlet.momentum.velocity.set_state(0.4)
18:     >>> inlet_turbulence = cold_inlet.turbulence
19:     >>> turbulence_specification = inlet_turbulence.turbulence_specification

## doc/source/user_guide/solver_settings/set_up/materials.rst
12:     >>> file_name = examples.download_file("mixing_elbow.cas.h5", "pyfluent/mixing_elbow")
13:     >>> solver_session = pyfluent.Solver.from_install()
14:     >>> solver_session.settings.file.read_case(file_name=file_name)
15:     >>> materials = pyfluent.Materials(settings_source=solver_session)
16:     >>> fluids = materials.fluid
17:     >>> fluids.make_a_copy(from_="air",to="air-2")
18:     >>> air_copy = fluids["air-2"]

## doc/source/user_guide/solver_settings/set_up/models/energy.rst
11:     >>> file_name = examples.download_file("mixing_elbow.cas.h5", "pyfluent/mixing_elbow")
12:     >>> solver_session = pyfluent.Solver.from_install()
13:     >>> solver_session.settings.file.read_case(file_name=file_name)
14:     >>> energy = pyfluent.solver.Energy(settings_source=solver_session)
15:     >>> energy.enabled.get_state()
16:     True
17:     >>> from pprint import pprint

## doc/source/user_guide/solver_settings/set_up/models/models_contents.rst
16:     >>> file_name = examples.download_file("mixing_elbow.cas.h5", "pyfluent/mixing_elbow")
17:     >>> solver_session = pyfluent.Solver.from_install()
18:     >>> solver_session.settings.file.read_case(file_name=file_name)
19:
20:
21: .. toctree::
22:    :maxdepth: 1

## doc/source/user_guide/solver_settings/solver_settings_contents.rst
366:   >>> import_file_name = examples.download_file("mixing_elbow.msh.h5", "pyfluent/mixing_elbow")
367:   >>> solver_session = pyfluent.Solver.from_install()
368:   >>> solver_session.settings.file.read(file_type="case", file_name=import_file_name)
369:   Fast-loading...
370:   ...Done
371:   >>> viscous = Viscous(settings_source=solver_session)
372:   >>> viscous.is_active()

## doc/source/user_guide/transfer_data.rst
18:   >>> mesh_file_name = download_file("mixing_elbow.msh.h5", "pyfluent/mixing_elbow")
19:
20:   >>> source = pyfluent.PureMeshing.from_install()
21:   >>> source.meshing_workflow.application.file.read_mesh(file_name=mesh_file_name)
22:
23:   >>> solver = pyfluent.Solver.from_install()
24:   >>> transfer_case(source_instance=source, solvers=[solver], file_type="mesh")
37:   >>> case_file_name = download_file("mixing_elbow.cas.h5", "pyfluent/mixing_elbow")
38:
39:   >>> source = pyfluent.Solver.from_install()
40:   >>> source.settings.file.read_case(file_name=case_file_name)
41:
42:   >>> target = pyfluent.Solver.from_install()
43:   >>> transfer_case(source_instance=source, solvers=[target])

## doc/source/user_guide/user_guide_contents.rst
45:   >>> watertight.import_geometry.file_name = pyfluent.examples.download_file("mixing_elbow.pmdb","pyfluent/mixing_elbow")
46:   >>> watertight.import_geometry()
47:   >>> watertight.create_volume_mesh()
48:   >>> solver_session = meshing_session.switch_to_solver()
49:   >>> setup, solution = solver_session.settings.setup, solver_session.settings.solution
50:   >>> setup.boundary_conditions.set_zone_type(zone_list=["cold-inlet", "hot-inlet"], new_type="velocity-inlet")
51:   >>> setup.boundary_conditions.set_zone_type(zone_list=["outlet"], new_type="pressure-outlet")

## examples/00-fluent/DOE_ML.py
85: import_filename = examples.download_file(
86:     "elbow.cas.h5",
87:     "pyfluent/examples/DOE-ML-Mixing-Elbow",
88:     save_path=os.getcwd(),
89: )
90:
91: #######################

## examples/00-fluent/Electrolysis_Modeling_workflow.py
106: mesh_file = examples.download_file(
107:     "electrolysis.msh.h5",
108:     "pyfluent/electrolysis",
109:     save_path=os.getcwd(),
110: )
111:
112: solver.settings.file.read_mesh(file_name=mesh_file)

## examples/00-fluent/Modeling_solidification_workflow.py
114: mesh_file = examples.download_file(
115:     "solid.msh",
116:     "pyfluent/solidification",
117:     save_path=os.getcwd(),
118: )
119:
120: solver.settings.file.read_mesh(file_name=mesh_file)

## examples/00-fluent/ahmed_body_workflow.py
91: geometry_filename = examples.download_file(
92:     filenames.get(platform.system(), filenames["Other"]),
93:     "pyfluent/examples/Ahmed-Body-Simulation",
94:     save_path=os.getcwd(),
95: )
96:
97: workflow.InitializeWorkflow(WorkflowType="Watertight Geometry")

## examples/00-fluent/axial_fan_performance_curve_workflow.py
52: import_file_name = examples.download_file(
53:     "axial_fan.msh.h5",
54:     "pyfluent/axial_fan_performance_curve",
55:     save_path=os.getcwd(),
56: )
57:
58: # %%

## examples/00-fluent/battery_pack.py
111: mesh_file = examples.download_file(
112:     "1P3S_battery_pack.msh",
113:     "pyfluent/battery_pack",
114:     save_path=os.getcwd(),
115: )
116:
117: solver.settings.file.read_mesh(file_name=mesh_file)

## examples/00-fluent/brake.py
61: import_filename = examples.download_file(
62:     "brake.msh.h5",
63:     "pyfluent/examples/Brake-Thermal-PyVista-Matplotlib",
64:     save_path=os.getcwd(),
65: )
66:
67: ####################################################################################

## examples/00-fluent/catalytic_converter_workflow.py
152: geometry_filename = examples.download_file(
153:     filenames.get(platform.system(), filenames["Other"]),
154:     "/pyfluent/catalytic_converter/",
155:     save_path=os.getcwd(),
156: )
157:
158: workflow.InitializeWorkflow(WorkflowType="Watertight Geometry")

## examples/00-fluent/conjugate_heat_transfer.py
81: geom_filename = examples.download_file(
82:     filenames.get(platform.system(), filenames["Other"]),
83:     "pyfluent/examples/CHT",
84:     save_path=os.getcwd(),
85: )
86:
87: #######################

## examples/00-fluent/exhaust_system_settings_api.py
74: import_file_name = examples.download_file(
75:     "exhaust_system.fmd", "pyfluent/exhaust_system"
76: )
77:
78: ###############################################################################
79: # Launch Fluent
80: # ~~~~~~~~~~~~~

## examples/00-fluent/external_compressible_flow.py
78:     examples.download_file(CAD_file, "pyfluent/external_compressible")
79:     for CAD_file in ["wing.scdoc", "wing.pmdb"]
80: ]
81:
82: ###############################################################################
83: # Launch Fluent
84: # ~~~~~~~~~~~~~

## examples/00-fluent/frozen_rotor_workflow.py
94: impeller_mesh = examples.download_file(
95:     "impeller.msh.h5",
96:     "pyfluent/examples/pump-volute",
97:     save_path=os.getcwd(),
98: )
99:
100: volute_mesh = examples.download_file(
101:     "volute.msh.h5",
102:     "pyfluent/examples/pump-volute",
103:     save_path=os.getcwd(),
104: )
105:
106: ################################################################################################################

## examples/00-fluent/fsi_1way_workflow.py
99: mesh_file = examples.download_file(
100:     "fsi_1way.msh.h5",
101:     "pyfluent/fsi_1way",
102:     save_path=os.getcwd(),
103: )
104: solver.settings.file.read_case(file_name=mesh_file)
105:

## examples/00-fluent/lunar_lander_thermal.py
105:     examples.download_file(
106:         f,
107:         "pyfluent/lunar_lander_thermal",
108:         save_path=os.getcwd(),
109:     )
110:     for f in [
111:         "lander_geom.scdoc",

## examples/00-fluent/mixing_elbow_settings_api.py
57: import_file_name = examples.download_file(
58:     "mixing_elbow.msh.h5",
59:     "pyfluent/mixing_elbow",
60:     save_path=os.getcwd(),
61: )
62:
63: ###############################################################################

## examples/00-fluent/mixing_tank_workflow.py
98: geometry_filename = examples.download_file(
99:     filenames.get(platform.system(), filenames["Other"]),
100:     "pyfluent/examples/MixingTank-WorkFlow",
101:     save_path=os.getcwd(),
102: )
103:
104: workflow.InitializeWorkflow(WorkflowType="Watertight Geometry")

## examples/00-fluent/modeling_ablation.py
93: import_filename = examples.download_file(
94:     "ablation.msh.h5",
95:     "pyfluent/examples/Ablation-tutorial",
96:     save_path=os.getcwd(),
97: )
98:
99: ####################################################################################

## examples/00-fluent/modeling_cavitation.py
73: cav_file = examples.download_file(
74:     "cav.msh.gz", "pyfluent/cavitation", save_path=os.getcwd()
75: )
76:
77: ###############################################################################
78: # Launch a Fluent session in the 2d solution mode with double precision running
79: # on four processors and print Fluent version.

## examples/00-fluent/parametric_static_mixer_1.py
59: import_filename = examples.download_file(
60:     "Static_Mixer_main.cas.h5",
61:     "pyfluent/static_mixer",
62:     save_path=os.getcwd(),
63: )
64:
65: solver_session.settings.file.read_case(file_name=import_filename)

## examples/00-fluent/radiation_headlamp.py
75:     examples.download_file(
76:         f,
77:         "pyfluent/radiation_headlamp",
78:         save_path=os.getcwd(),
79:     )
80:     for f in ["headlamp.scdoc", "headlamp.pmdb"]
81: ]

## examples/00-fluent/single_battery_cell_workflow.py
63: unit_battery_mesh = examples.download_file(
64:     "unit_battery.msh.h5",
65:     "pyfluent/battery_thermal_simulation",
66:     save_path=os.getcwd(),
67: )
68:
69: #######################################################################################

## examples/00-fluent/species_transport.py
118: mesh_file = download_file(
119:     "gascomb.msh.gz", "pyfluent/tutorials/species_transport", save_path=os.getcwd()
120: )
121: solver_session.settings.file.read_mesh(file_name=mesh_file)
122:
123: # %%
124: # General Settings

## examples/00-fluent/steady_vortex.py
131: vortex_mesh = download_file(
132:     "vortex-mixingtank.msh.h5",
133:     "pyfluent/examples/Steady-Vortex-VOF",
134:     save_path=os.getcwd(),
135: )
136:
137: # %%

## examples/00-fluent/transient_compressible_nozzle_workflow.py
111: geometry_filename = examples.download_file(
112:     filenames.get(platform.system(), filenames["Other"]),
113:     "pyfluent/transient_compressible_simulation",
114:     save_path=os.getcwd(),
115: )
116:
117: workflow.InitializeWorkflow(WorkflowType="Watertight Geometry")

## examples/00-fluent/tyler_sofrin_modes.py
130: import_filename = examples.download_file(
131:     "axial_comp_fullWheel_DFT_23R2.cas.h5",
132:     "pyfluent/examples/Tyler-Sofrin-Modes-Compressor",
133:     save_path=os.getcwd(),
134: )
135:
136: examples.download_file(
137:     "axial_comp_fullWheel_DFT_23R2.dat.h5",
138:     "pyfluent/examples/Tyler-Sofrin-Modes-Compressor",
139:     save_path=os.getcwd(),
140: )
141:
142: #######################################################################################

## src/ansys/fluent/core/examples/downloads.py
121: def download_file(
122:     file_name: str,
123:     directory: str | None = None,
124:     save_path: "PathType | None" = None,
125:     return_without_path: bool | None = None,
126: ) -> str:
127:     """Download specified example file from the Ansys example data repository.
160:     >>> file_path = examples.download_file("bracket.iges", "geometry")
161:     >>> file_path
162:     '/home/user/.local/share/ansys_fluent_core/examples/bracket.iges'
163:     >>> file_name = examples.download_file("bracket.iges", "geometry", return_without_path=True)
164:     >>> file_name
165:     'bracket.iges'
166:     >>> file_path = examples.download_file("bracket.iges", "geometry", save_path='.')
167:     '/home/<current_folder_path>/bracket.iges'
168:     >>> file_name = examples.download_file("bracket.iges", "geometry", save_path='.', return_without_path=True)
169:     >>> file_name
170:     'bracket.iges'
171:     >>> file_path = examples.download_file("bracket.iges", "geometry", save_path='<user_specified_path>')
172:     '/home/<user_specified_path>/bracket.iges'
173:     >>> file_name = examples.download_file("bracket.iges", "geometry", save_path='<user_specified_path>',
174:     ...                                   return_without_path=True)
175:     >>> file_name
176:     'bracket.iges'
177:     """
178:     if return_without_path is None:
179:         if pyfluent.config.launch_fluent_container:
