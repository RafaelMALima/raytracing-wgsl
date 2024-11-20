fn hit_sphere(center: vec3f, radius: f32, r: ray, record: ptr<function, hit_record>, max: f32)
{
  var t: f32 = 0.;
  var p0 : vec3f = r.origin;
  var d : vec3f = r.direction;
  //calcula o baskhara para a interseccao
  var a = dot(d,d);
  var b = 2*dot(d,p0-center);
  var c = dot(p0 - center,p0 - center) - radius*radius;
  var delta  = b*b - 4*a*c;
  if (delta > 0){
    record.hit_anything = false;
  }
  if (delta > 0){
    var sqrt_delta = sqrt(delta);
    var t1: f32 = (-b - sqrt_delta) / (2.0 * a);
    var t2: f32 = (-b + sqrt_delta) / (2.0 * a);

    // Initialize t to maximum value
    var t_min = 0.01;  // Small epsilon to avoid self-intersection
    var t = max;       // Use the 'max' parameter provided

    var invert_normal = 1.;
    // Find the nearest positive t
    if (t1 > t_min && t1 < t) {
      t = t1;
      record.frontface = true;
    }
    if (t2 > t_min && t2 < t) {
      t = t2;
      record.frontface = false;
      invert_normal = -1.;
    }

    // If t remains unchanged, no valid intersection was found
    if (t == max) {
      record.hit_anything = false;
      return;
    }

    // Record the hit information
    record.hit_anything = true;
    record.t = t;
    record.p = p0 + t * d;

    record.normal = normalize(record.p - center)*invert_normal;
  }
}

fn hit_quad(r: ray, Q: vec4f, u: vec4f, v: vec4f, record: ptr<function, hit_record>, max: f32)
{
  var n = cross(u.xyz, v.xyz);
  var normal = normalize(n);
  var D = dot(normal, Q.xyz);
  var w = n / dot(n.xyz, n.xyz);

  var denom = dot(normal, r.direction);
  if (abs(denom) < 0.0001)
  {
    record.hit_anything = false;
    return;
  }

  var t = (D - dot(normal, r.origin)) / denom;
  if (t < RAY_TMIN || t > max)
  {
    record.hit_anything = false;
    return;
  }

  var intersection = ray_at(r, t);
  var planar_hitpt_vector = intersection - Q.xyz;
  var alpha = dot(w, cross(planar_hitpt_vector, v.xyz));
  var beta = dot(w, cross(u.xyz, planar_hitpt_vector));

  if (alpha < 0.0 || alpha > 1.0 || beta < 0.0 || beta > 1.0)
  {
    record.hit_anything = false;
    return;
  }

  if (dot(normal, r.direction) > 0.0)
  {
    record.hit_anything = false;
    return;
  }

  record.t = t;
  record.p = intersection;
  record.normal = normal;
  record.hit_anything = true;
}

fn hit_triangle(r: ray, v0: vec3f, v1: vec3f, v2: vec3f, record: ptr<function, hit_record>, max: f32)
{
  var v1v0 = v1 - v0;
  var v2v0 = v2 - v0;
  var rov0 = r.origin - v0;

  var n = cross(v1v0, v2v0);
  var q = cross(rov0, r.direction);

  var d = 1.0 / dot(r.direction, n);

  var u = d * dot(-q, v2v0);
  var v = d * dot(q, v1v0);
  var t = d * dot(-n, rov0);

  if (u < 0.0 || u > 1.0 || v < 0.0 || (u + v) > 1.0)
  {
    record.hit_anything = false;
    return;
  }

  if (t < RAY_TMIN || t > max)
  {
    record.hit_anything = false;
    return;
  }

  record.t = t;
  record.p = ray_at(r, t);
  record.normal = normalize(n);
  record.hit_anything = true;
}

fn hit_box(r: ray, center: vec3f, rad: vec3f, record: ptr<function, hit_record>, t_max: f32)
{
  var m = 1.0 / r.direction;
  var n = m * (r.origin - center);
  var k = abs(m) * rad;

  var t1 = -n - k;
  var t2 = -n + k;

  var tN = max(max(t1.x, t1.y), t1.z);
  var tF = min(min(t2.x, t2.y), t2.z);

  if (tN > tF || tF < 0.0)
  {
    record.hit_anything = false;
    return;
  }

  var t = tN;
  if (t < RAY_TMIN || t > t_max)
  {
    record.hit_anything = false;
    return;
  }

  record.t = t;
  record.p = ray_at(r, t);
  record.normal = -sign(r.direction) * step(t1.yzx, t1.xyz) * step(t1.zxy, t1.xyz);
  record.hit_anything = true;

  return;
}

fn hit_plane(r: ray, n: vec3f, offset: f32, record: ptr<function, hit_record>, max: f32)
{
  var normal = normalize(n);
  if abs(dot(normal, r.origin)) < 0.001 {
    record.hit_anything=false;
    return;
  }
  var t = -(dot(normal, r.origin) + offset) / dot(normal, r.direction);
    if (t < RAY_TMIN || t > max)
  {
    record.hit_anything=false;
    return;
  }    
  record.t = t;
  record.p = r.origin + t * r.direction;
  record.normal = normal;
  record.hit_anything = true;
}