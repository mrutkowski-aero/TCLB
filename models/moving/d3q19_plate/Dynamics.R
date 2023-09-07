MRTMAT = matrix(c(
1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
-30,-11,-11,-11,-11,-11,-11,8,8,8,8,8,8,8,8,8,8,8,8,
12,-4,-4,-4,-4,-4,-4,1,1,1,1,1,1,1,1,1,1,1,1,
0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1,0,0,0,0,
0,-4,4,0,0,0,0,1,-1,1,-1,1,-1,1,-1,0,0,0,0,
0,0,0,1,-1,0,0,1,1,-1,-1,0,0,0,0,1,-1,1,-1,
0,0,0,-4,4,0,0,1,1,-1,-1,0,0,0,0,1,-1,1,-1,
0,0,0,0,0,1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1,
0,0,0,0,0,-4,4,0,0,0,0,1,1,-1,-1,1,1,-1,-1,
 0,2,2,-1,-1,-1,-1,1,1,1,1,1,1,1,1,-2,-2,-2,-2,
0,-4,-4,2,2,2,2,1,1,1,1,1,1,1,1,-2,-2,-2,-2,
 0,0,0,1,1,-1,-1,1,1,1,1,-1,-1,-1,-1,0,0,0,0,
0,0,0,-2,-2,2,2,1,1,1,1,-1,-1,-1,-1,0,0,0,0,
 0,0,0,0,0,0,0,1,-1,-1,1,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,-1,1,
 0,0,0,0,0,0,0,0,0,0,0,1,-1,-1,1,0,0,0,0,
0,0,0,0,0,0,0,1,-1,1,-1,-1,1,-1,1,0,0,0,0,
0,0,0,0,0,0,0,-1,-1,1,1,0,0,0,0,1,-1,1,-1,
0,0,0,0,0,0,0,0,0,0,0,1,1,-1,-1,-1,-1,1,1
),19,19)

v = diag(t(MRTMAT) %*% MRTMAT)
MRTMAT.inv = diag(1/v) %*% t(MRTMAT)

selU = c(4,6,8)

AddDensity(
	name = paste("f",0:18,sep=""),
	dx   = MRTMAT[,selU[1]],
	dy   = MRTMAT[,selU[2]],
	dz   = MRTMAT[,selU[3]],
	comment=paste("density F",0:18),
	group="f"
)

# Boundary initialization

#AddDensity( name="avg_ux", group="avg_u")
#AddDensity( name="avg_uy", group="avg_u")
#AddDensity( name="avg_uz", group="avg_u")
#AddDensity( name="avg_fx", group="avg_f")
#AddDensity( name="avg_fy", group="avg_f")
#AddDensity( name="avg_fz", group="avg_f")

AddNodeType(name="WVelocity", group="BOUNDARY")
AddNodeType(name="EPressure", group="BOUNDARY")
AddNodeType(name="Solid", group="BOUNDARY")
AddNodeType(name="Wall", group="BOUNDARY")
AddNodeType(name="MRT", group="COLLISION")

# Quantities - table of fields that can be exported from the LB lattice (like density, velocity etc)
#  name - name of the field
#  type - C type of the field, "real_t" - for single/double float, and "vector_t" for 3D vector single/double float
# Every field must correspond to a function in "Dynamics.c".
# If one have filed [something] with type [type], one have to define a function: 
# [type] get[something]() { return ...; }

AddQuantity(name="Rho",unit="kg/m3")
AddQuantity(name="U",unit="m/s",vector=T)
#AddQuantity(name="U_AVG",unit="m/s",vector=T)
#AddQuantity(name="F_AVG",unit="N/m3",vector=T)
#AddQuantity(name="Solid",unit="1")
AddQuantity( name="RhoB",adjoint=T)
AddQuantity( name="UB",adjoint=T,vector=T)

# Settings - table of settings (constants) that are taken from a .xml file
#  name - name of the constant variable
#  comment - additional comment
# You can state that another setting is 'derived' from this one stating for example: omega='1.0/(3*nu + 0.5)'

AddSetting(name="omega", comment='one over relaxation time')
AddSetting(name="nu", omega='1.0/(3*nu + 0.5)', default=0.16666666, comment='viscosity')
#AddSetting(name="nu", default=0.16666666, comment='viscosity', zonal=T, unit="m2/s")
AddSetting(name="Velocity", default=0, comment='inlet/outlet/init velocity', zonal=T, unit="m/s")
#AddSetting(name="VelocityX", default=0, comment='inlet/outlet/init velocity', zonal=T, unit="m/s")
#AddSetting(name="VelocityY", default=0, comment='inlet/outlet/init velocity', zonal=T, unit="m/s")
#AddSetting(name="VelocityZ", default=0, comment='inlet/outlet/init velocity', zonal=T, unit="m/s")
AddSetting(name="Density", default=1, comment='inlet/outlet/init density', zonal=T, unit="kg/m3")
AddSetting(name="Smag", default=1, comment='inlet density')
#AddSetting(name="ExternalForceX", default=0, comment='external force x', zonal=T, unit="N/m3")
#AddSetting(name="ExternalForceY", default=0, comment='external force y', zonal=T, unit="N/m3")
#AddSetting(name="ExternalForceZ", default=0, comment='external force z', zonal=T, unit="N/m3")
#AddSetting(name="PDX", default=0, comment='plate dimension X', unit="m")
#AddSetting(name="PDY", default=0, comment='plate dimension Y', unit="m")
#AddSetting(name="PDZ", default=0, comment='plate dimension Z', unit="m")
#AddSetting(name="PRAD", default=0, comment='cylinder radious', unit='m')
#AddSetting(name="SM",   default=1, comment='smoothing diameter', unit="m")
#AddSetting(name="SM_M", default=0, comment='smoothing bias')
#AddSetting(name="EPSF", default=1, comment='boundary function, 0 - linear boundary,1 - #third order boundary')
#AddSetting(name="BF", default=0, comment='beta function bool')
#AddSetting(name="PX", default=0, comment='plate position X', zonal=T, unit="m")
#AddSetting(name="PY", default=0, comment='plate position Y', zonal=T, unit="m")
#AddSetting(name="PZ", default=0, comment='plate position Z', zonal=T, unit="m")
#AddSetting(name="PRX", default=0, comment='plate angle', zonal=T, unit="1")
#AddSetting(name="PRY", default=0, comment='plate angle', zonal=T, unit="1")
#AddSetting(name="PRZ", default=0, comment='plate angle', zonal=T, unit="1")


# Globals - table of global integrals that can be monitored and optimized

#AddGlobal(name="ForceX", comment='reaction force X', unit="N/m")
#AddGlobal(name="ForceY", comment='reaction force Y', unit="N/m")
#AddGlobal(name="ForceZ", comment='reaction force Y', unit="N/m")
#AddGlobal(name="Moment", comment='reaction moment', unit="N")
#AddGlobal(name="PowerX", comment='Translational power X', unit="W/m")
#AddGlobal(name="PowerY", comment='Translational power Y', unit="W/m")
#AddGlobal(name="PowerZ", comment='Translational power Z', unit="W/m")
#AddGlobal(name="PowerR", comment='Rotational Power', unit="W/m")
#AddGlobal(name="Power", comment='Fluids power', unit="W/m")
#AddGlobal(name="Power2", comment='Fluids power 2', unit="W/m")
#AddGlobal(name="VolumeW", comment="Volume of moving body", unit="m2")




#AddObjective("EfficiencyX", PV("ForceX") * PV("Power") ^ (-1))
#AddObjective("EfficiencyY", PV("ForceY") * PV("Power") ^ (-1))
#AddObjective("EfficiencyZ", PV("ForceZ") * PV("Power") ^ (-1))



