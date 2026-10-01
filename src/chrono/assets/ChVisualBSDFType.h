#ifndef CHVISUALBSDFTYPE_H
#define CHVISUALBSDFTYPE_H

enum class BSDFType {
	DIFFUSE,
	SPECULAR,
	DIELECTRIC,
	GLOSSY,
	SIMPLEPRINCIPLED,
	PRINCIPLED,
	HAPKE,
	RETROREFLECTIVE,
	VDB,
	VDBHAPKE,
	VDBVOL,
	PLANET  // a planet seen from space (Vulkan RT): its day map, night lights, water glint and clouds, under the scene's atmosphere
};

#endif // CHVISUALBSDFTYPE_H