const THREAD_COUNT = 16;
const RAY_TMIN = 0.0001;
const RAY_TMAX = 100.0;
const PI = 3.1415927f;
const FRAC_1_PI = 0.31830987f;
const FRAC_2_PI = 1.5707964f;

@group(0) @binding(0)  
  var<storage, read_write> fb : array<vec4f>;

@group(0) @binding(1)
  var<storage, read_write> rtfb : array<vec4f>;

@group(1) @binding(0)
  var<storage, read_write> uniforms : array<f32>;

@group(2) @binding(0)
  var<storage, read_write> spheresb : array<sphere>;

@group(2) @binding(1)
  var<storage, read_write> quadsb : array<quad>;

@group(2) @binding(2)
  var<storage, read_write> boxesb : array<box>;

@group(2) @binding(3)
  var<storage, read_write> trianglesb : array<triangle>;

@group(2) @binding(4)
  var<storage, read_write> meshb : array<mesh>;

//@group(2) @binding(5)
//  var<storage, read_write> planesb : array<plane>;

struct ray {
  origin : vec3f,
  direction : vec3f,
};

struct sphere {
  transform : vec4f,
  color : vec4f,
  material : vec4f,
};

struct quad {
  Q : vec4f,
  u : vec4f,
  v : vec4f,
  color : vec4f,
  material : vec4f,
};

struct box {
  center : vec4f,
  radius : vec4f,
  rotation: vec4f,
  color : vec4f,
  material : vec4f,
};

struct triangle {
  v0 : vec4f,
  v1 : vec4f,
  v2 : vec4f,
};

struct mesh {
  transform : vec4f,
  scale : vec4f,
  rotation : vec4f,
  color : vec4f,
  material : vec4f,
  min : vec4f,
  max : vec4f,
  show_bb : f32,
  start : f32,
  end : f32,
};

struct plane {
  normal : vec4f,
  color : vec4f,
  material : vec4f,
  offset : f32,
};


struct material_behaviour {
  scatter : bool,
  direction : vec3f,
};

struct camera {
  origin : vec3f,
  lower_left_corner : vec3f,
  horizontal : vec3f,
  vertical : vec3f,
  u : vec3f,
  v : vec3f,
  w : vec3f,
  lens_radius : f32,
};

struct hit_record {
  t : f32,
  p : vec3f,
  normal : vec3f,
  object_color : vec4f,
  object_material : vec4f,
  frontface : bool,
  hit_anything : bool,
};

fn ray_at(r: ray, t: f32) -> vec3f
{
  return r.origin + t * r.direction;
}

fn get_ray(cam: camera, uv: vec2f, rng_state: ptr<function, u32>) -> ray
{
  var rd = cam.lens_radius * rng_next_vec3_in_unit_disk(rng_state);
  var offset = cam.u * rd.x + cam.v * rd.y;
  return ray(cam.origin + offset, normalize(cam.lower_left_corner + uv.x * cam.horizontal + uv.y * cam.vertical - cam.origin - offset));
}

fn get_camera(lookfrom: vec3f, lookat: vec3f, vup: vec3f, vfov: f32, aspect_ratio: f32, aperture: f32, focus_dist: f32) -> camera
{
  var camera = camera();
  camera.lens_radius = aperture / 2.0;

  var theta = degrees_to_radians(vfov);
  var h = tan(theta / 2.0);
  var w = aspect_ratio * h;

  camera.origin = lookfrom;
  camera.w = normalize(lookfrom - lookat);
  camera.u = normalize(cross(vup, camera.w));
  camera.v = cross(camera.u, camera.w);

  camera.lower_left_corner = camera.origin - w * focus_dist * camera.u - h * focus_dist * camera.v - focus_dist * camera.w;
  camera.horizontal = 2.0 * w * focus_dist * camera.u;
  camera.vertical = 2.0 * h * focus_dist * camera.v;

  return camera;
}

fn environment_color(direction: vec3f, color1: vec3f, color2: vec3f) -> vec3f
{
  var unit_direction = normalize(direction);
  var t = 0.5 * (unit_direction.y + 1.0);
  var col = (1.0 - t) * color1 + t * color2;

  var sun_direction = normalize(vec3(uniforms[13], uniforms[14], uniforms[15]));
  var sun_color = int_to_rgb(i32(uniforms[17]));
  var sun_intensity = uniforms[16];
  var sun_size = uniforms[18];

  var sun = clamp(dot(sun_direction, unit_direction), 0.0, 1.0);
  col += sun_color * max(0, (pow(sun, sun_size) * sun_intensity));

  return col;
}

fn check_ray_collision(r: ray, max: f32) -> hit_record
{
  var spheresCount = i32(uniforms[19]);
  var quadsCount = i32(uniforms[20]);
  var boxesCount = i32(uniforms[21]);
  var trianglesCount = i32(uniforms[22]);
  var meshCount = i32(uniforms[27]);
  //var planeCount = i32(uniforms[28]);

  var record = hit_record(RAY_TMAX, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);
  var closest = record;
  for (var i = 0; i < spheresCount; i++){
    var sphere = spheresb[i];
    var record = hit_record(RAY_TMAX, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);
    hit_sphere(sphere.transform.xyz, sphere.transform.w, r, &record, max);
    //check if this record is closest hit if hit happens
    if (record.hit_anything && (closest.hit_anything == false || length(closest.p - r.origin) > length(record.p - r.origin))){
      record.object_color = sphere.color;
      record.object_material = sphere.material;
      closest = record;
    }
  }
  for (var i = 0; i < boxesCount; i++){
    var box = boxesb[i];
    var record = hit_record(RAY_TMAX, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);
    //ja que nao tem mais espaco para mais um buffer, gambiarra aqui para a primitiva do plano funcionar

    //se tem a flag, eh um plano
    if (box.center.w > 1){
      var plane = box;
      hit_plane(r, plane.center.xyz, plane.radius.x, &record, max);
      if (record.hit_anything){
        record.object_color = plane.color;
        record.object_material = plane.material;
        closest = record;
      }
    }
    //se nao tem a flag, eh uma box normal
    else{
      hit_box(r, box.center.xyz, box.radius.xyz, &record, max);
      //check if this record is closest hit if hit happens
      if (record.hit_anything && (closest.hit_anything == false || length(closest.p - r.origin) > length(record.p - r.origin))){
        record.object_color = box.color;
        record.object_material = box.material;
        closest = record;
      }
    }
  }
  for (var i = 0; i < quadsCount; i++){
    var quad = quadsb[i];
    var record = hit_record(RAY_TMAX, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);
    hit_quad(r, quad.Q, quad.u, quad.v, &record, max);
    if (record.hit_anything && (closest.hit_anything == false || length(closest.p - r.origin) > length(record.p - r.origin))){
      record.object_color = quad.color;
      record.object_material = quad.material;
      closest = record;
    }
  }
  for (var i = 0; i < meshCount; i++){
    var mesh = meshb[i];
    for (var j = i32(mesh.start); j < i32(mesh.end); j++){
      var tri = trianglesb[j];
      var mesh_position = mesh.transform.xyz;
      hit_triangle(r, mesh_position + tri.v0.xyz,mesh_position+  tri.v1.xyz,mesh_position+ tri.v2.xyz,&record, max);
      if (record.hit_anything && (closest.hit_anything == false || length(closest.p - r.origin) > length(record.p - r.origin))){
        record.object_color = mesh.color;
        record.object_material = mesh.material;
        closest = record;
      }

    }
  }

  return closest;
}

fn lambertian(normal : vec3f, absorption: f32, random_sphere: vec3f, rng_state: ptr<function, u32>) -> material_behaviour
{
  var dir = random_sphere + normal;
  return material_behaviour(true, normalize(dir));
}

fn metal(normal : vec3f, direction: vec3f, fuzz: f32, random_sphere: vec3f) -> material_behaviour
{
  var reflect = direction - 2* dot(normal,direction)*normal;
  var lamb_dir = random_sphere + normal;
  var reflect_direction = normalize(reflect)+fuzz*lamb_dir;
  return material_behaviour(true, normalize(reflect_direction));
}

fn schlick(cosine: f32, refraction_index: f32) -> f32
{
  var r0 = (1.0 - refraction_index) / (1.0 + refraction_index);
  r0 = r0 * r0;
  return r0 + (1.0 - r0) * pow(1.0 - cosine, 5.0);
}

fn schlick_aproximation(refraction_index: f32, cos_theta: f32) -> f32{
  var r0 = pow((1 - refraction_index)/(1 + refraction_index), 2.);
  var r_theta = r0 + (1-r0)*pow((1-cos_theta), 5.);
  return r_theta;
}

fn dielectric(normal : vec3f, r_direction: vec3f, refraction_index: f32, frontface: bool, random_sphere: vec3f, fuzz: f32, rng_state: ptr<function, u32>) -> material_behaviour{    
  let unit_direction = normalize(r_direction);
  let cos_theta = min(dot(-unit_direction, normal), 1.0);
  let sin_theta = sqrt(1.0 - cos_theta * cos_theta);


  var refraction_ratio = 0.;
  if frontface { refraction_ratio = 1.0 / refraction_index; } else { refraction_ratio = refraction_index; };

  var r_theta = schlick_aproximation(refraction_ratio, refraction_ratio);
  var direction: vec3f;
  if (refraction_ratio * sin_theta > 1.0 || r_theta > rng_next_float(rng_state)) {
    // Total internal reflection
    direction = reflect(unit_direction, normal);
  } else {
      // Refraction
      direction = refract(unit_direction, normal, refraction_ratio);
  }

  // Apply fuzziness if needed
  return material_behaviour(true, direction + fuzz * random_sphere);
}


fn emmisive(color: vec3f, light: f32) -> material_behaviour
{
  return material_behaviour(false, vec3f(0.0));
}

fn trace(r: ray, rng_state: ptr<function, u32>) -> vec3f
{
  var maxbounces = i32(uniforms[2]);
  var light = vec3f(0.0);
  var color = vec3f(1.0);
  var r_ = r;
  
  var backgroundcolor1 = int_to_rgb(i32(uniforms[11]));
  var backgroundcolor2 = int_to_rgb(i32(uniforms[12]));
  var behaviour = material_behaviour(true, vec3f(0.0));

  var accumulated_color = color;
  var hitrec = hit_record(RAY_TMAX, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);
  var spec_effect = 1.;

  for (var j = 0; j < maxbounces; j = j + 1)
  {
    hitrec = check_ray_collision(r_, 1000.);
    if (!hitrec.hit_anything){ light = accumulated_color *environment_color(r_.direction, backgroundcolor1, backgroundcolor2); break; }
    //determine bounce behaviour
    var smoothn = hitrec.object_material.x;
    var absorp = hitrec.object_material.y;
    var spec = hitrec.object_material.z;
    var luce = hitrec.object_material.w;

    if (luce > 0){
      light = accumulated_color * luce*hitrec.object_color.xyz;
      break;
    }
    if (smoothn > 0){
      //deal with specular
      if (rng_next_float(rng_state) > spec){
        behaviour = lambertian(hitrec.normal, absorp, rng_next_vec3_in_unit_sphere(rng_state), rng_state);
        accumulated_color *= hitrec.object_color.xyz*(1-absorp);
      } else  {
        behaviour = metal(hitrec.normal, r_.direction, 1-smoothn*(1-absorp), rng_next_vec3_in_unit_sphere(rng_state));
        //var new_color = hitrec.object_color.xyz;
        //accumulated_color *= new_color + smoothn*spec*(vec3f(1) - new_color);
      }
    }
    else if (smoothn < 0){
      behaviour = dielectric(hitrec.normal, r_.direction, spec, hitrec.frontface, rng_next_vec3_in_unit_sphere(rng_state), absorp, rng_state);
      //behaviour = lambertian(hitrec.normal, absorp, rng_next_vec3_in_unit_sphere(rng_state), rng_state);
    }
    else {
      behaviour = lambertian(hitrec.normal, absorp, rng_next_vec3_in_unit_sphere(rng_state), rng_state);
      accumulated_color *= hitrec.object_color.xyz*(1-absorp);
    }
    if (behaviour.scatter){
      //accumulated_color *= hitrec.object_color.xyz;
      r_.origin = hitrec.p;
      r_.direction = behaviour.direction;
    }
    else{
      //accumulated_color *= hitrec.object_color.xyz;
      break;
    }
  }
  return light;
}

@compute @workgroup_size(THREAD_COUNT, THREAD_COUNT, 1)
fn render(@builtin(global_invocation_id) id : vec3u)
{
    var rez = uniforms[1];
    var time = u32(uniforms[0]);

    // init_rng (random number generator) we pass the pixel position, resolution and frame
    var rng_state = init_rng(vec2(id.x, id.y), vec2(u32(rez)), time);

    // Get uv
    var fragCoord = vec2f(f32(id.x), f32(id.y));
    var uv = (fragCoord + sample_square(&rng_state)) / vec2(rez);

    // Camera
    var lookfrom = vec3(uniforms[7], uniforms[8], uniforms[9]);
    var lookat = vec3(uniforms[23], uniforms[24], uniforms[25]);

    // Get camera
    var cam = get_camera(lookfrom, lookat, vec3(0.0, 1.0, 0.0), uniforms[10], 1.0, uniforms[6], uniforms[5]);
    var samples_per_pixel = i32(uniforms[4]);

    var color = vec3f(0f);

    // Steps:
    // 1. Loop for each sample per pixel
    for (var i = 0i; i < samples_per_pixel; i++){
      // 2. Get ray
      var point_pos = vec3f(uv - 0.5,1);
      var r = get_ray(cam, uv, &rng_state);
      // 3. Call trace function
      var ray_color = trace(r, &rng_state);
      color += ray_color;
    }
    // 4. Average the color
    color = color / f32(samples_per_pixel);

    var color_out = vec4(linear_to_gamma(color), 1.0);
    var map_fb = mapfb(id.xy, rez);
    
    // 5. Accumulate the color
    var should_accumulate = uniforms[3];

    var previous_color = rtfb[map_fb];
    var accumulated_color =  previous_color * should_accumulate + color_out;


    // Set the color to the framebuffer
    rtfb[map_fb] = accumulated_color;
    fb[map_fb] = accumulated_color / accumulated_color.w;
}