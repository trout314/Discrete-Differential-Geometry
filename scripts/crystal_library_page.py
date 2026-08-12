#!/usr/bin/env python3
"""Render the TCP crystal library as one self-contained HTML catalog.

Consumes data/crystal_gas/library.json (scripts/crystal_library_data.py) and
emits a single file with no external requests: the 3D viewers are a small
hand-written canvas renderer, the geometry is inlined as flat integer arrays.

Usage:
    python scripts/crystal_library_page.py --out /path/to/page.html
"""
import argparse
import html
import json
import math
import os

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ORDER = ["a15", "sigma", "z", "p", "delta", "r", "mu", "c15", "c14", "c36"]
NICE = {"a15": "A15", "c15": "C15", "c14": "C14", "c36": "C36",
        "sigma": "σ", "z": "Z", "mu": "μ", "r": "R", "p": "P",
        "delta": "δ"}
ZC = ["Z12", "Z14", "Z15", "Z16"]

CSS = """
*,*::before,*::after{box-sizing:border-box}
:root{
  --paper:#f2f4f6; --panel:#fbfcfd; --ink:#151a20; --dim:#5b6570;
  --faint:#8b95a0; --rule:#d7dde3; --rule2:#e7ecf0;
  --six:#c17d1e; --six-soft:#e8c98d; --five:#98a3ae;
  --z12:#9db1c4; --z14:#5b86a7; --z15:#356480; --z16:#1d4257;
  --ok:#2f6b4f; --warn:#a86b1c; --bad:#9c3b32;
  --onsix:#ffffff;
  --shadow:0 1px 2px rgba(20,30,40,.06),0 8px 24px -12px rgba(20,30,40,.18);
}
@media (prefers-color-scheme:dark){
  :root{
    --paper:#0e1216; --panel:#151b21; --ink:#dfe5ea; --dim:#94a0ac;
    --faint:#6e7a86; --rule:#242c34; --rule2:#1b2229;
    --six:#e0a54a; --six-soft:#7a5c26; --five:#4d5866;
    --z12:#7f95aa; --z14:#5b86a7; --z15:#78aecd; --z16:#a9cde3;
    --ok:#6fbf95; --warn:#d9a24e; --bad:#e08078; --onsix:#12171c;
    --shadow:0 1px 2px rgba(0,0,0,.4),0 10px 30px -14px rgba(0,0,0,.7);
  }
}
:root[data-theme="dark"]{
  --paper:#0e1216; --panel:#151b21; --ink:#dfe5ea; --dim:#94a0ac;
  --faint:#6e7a86; --rule:#242c34; --rule2:#1b2229;
  --six:#e0a54a; --six-soft:#7a5c26; --five:#4d5866;
  --z12:#7f95aa; --z14:#5b86a7; --z15:#78aecd; --z16:#a9cde3;
  --ok:#6fbf95; --warn:#d9a24e; --bad:#e08078; --onsix:#12171c;
  --shadow:0 1px 2px rgba(0,0,0,.4),0 10px 30px -14px rgba(0,0,0,.7);
}
:root[data-theme="light"]{
  --paper:#f2f4f6; --panel:#fbfcfd; --ink:#151a20; --dim:#5b6570;
  --faint:#8b95a0; --rule:#d7dde3; --rule2:#e7ecf0;
  --six:#c17d1e; --six-soft:#e8c98d; --five:#98a3ae;
  --z12:#9db1c4; --z14:#5b86a7; --z15:#356480; --z16:#1d4257;
  --ok:#2f6b4f; --warn:#a86b1c; --bad:#9c3b32; --onsix:#ffffff;
  --shadow:0 1px 2px rgba(20,30,40,.06),0 8px 24px -12px rgba(20,30,40,.18);
}
html{-webkit-text-size-adjust:100%}
body{margin:0;background:var(--paper);color:var(--ink);
  font:400 16px/1.6 ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,
       "Helvetica Neue",Arial,sans-serif;
  font-feature-settings:"kern" 1}
.serif{font-family:ui-serif,"Iowan Old Style","Palatino Linotype",Palatino,
  "Book Antiqua",Georgia,serif}
.mono,td.n,th.n{font-family:ui-monospace,SFMono-Regular,"SF Mono",Menlo,
  Consolas,"Liberation Mono",monospace;font-variant-numeric:tabular-nums}
.wrap{max-width:1180px;margin:0 auto;padding:0 28px}
a{color:inherit;text-decoration:none;border-bottom:1px solid var(--six-soft)}
a:hover{border-bottom-color:var(--six)}
:focus-visible{outline:2px solid var(--six);outline-offset:2px;border-radius:3px}

/* masthead */
header.mast{border-bottom:1px solid var(--rule);padding:56px 0 34px;
  background:linear-gradient(180deg,var(--panel),var(--paper))}
.eyebrow{font-size:11.5px;letter-spacing:.16em;text-transform:uppercase;
  color:var(--faint);margin:0 0 14px}
h1{font-size:clamp(30px,5vw,50px);line-height:1.06;margin:0 0 18px;
  font-weight:600;letter-spacing:-.015em;text-wrap:balance}
.lede{max-width:64ch;color:var(--dim);font-size:17px;margin:0 0 26px}
.thesis{display:flex;flex-wrap:wrap;gap:26px 44px;padding:20px 24px;
  border:1px solid var(--rule);border-left:3px solid var(--six);
  background:var(--panel);border-radius:2px}
.thesis dt{font-size:11.5px;letter-spacing:.13em;text-transform:uppercase;
  color:var(--faint);margin-bottom:5px}
.thesis dd{margin:0;font-size:19px;letter-spacing:-.01em}
.thesis .sm{font-size:13px;color:var(--dim);margin-top:3px;letter-spacing:0}

/* sections */
section{padding:46px 0;border-bottom:1px solid var(--rule2)}
h2{font-size:12px;letter-spacing:.17em;text-transform:uppercase;
  color:var(--faint);margin:0 0 6px;font-weight:600}
h2+.h2sub{font-size:23px;margin:0 0 22px;letter-spacing:-.012em;
  font-weight:600;text-wrap:balance}
p.note{max-width:68ch;color:var(--dim);font-size:14.5px}

/* tables */
.scroll{overflow-x:auto;-webkit-overflow-scrolling:touch}
table{border-collapse:collapse;width:100%;font-size:13.5px}
th,td{padding:7px 11px;text-align:right;white-space:nowrap;
  border-bottom:1px solid var(--rule2)}
th{font-size:10.5px;letter-spacing:.08em;text-transform:uppercase;
  color:var(--faint);font-weight:600;border-bottom:1px solid var(--rule);
  vertical-align:bottom}
th:first-child,td:first-child{text-align:left}
tbody tr:hover{background:var(--panel)}
td.name{font-weight:600}
td.name span.p{display:block;font-size:11.5px;color:var(--faint);
  font-weight:400}
.grp{border-left:1px solid var(--rule)}
tr.sep td{border-bottom:2px solid var(--rule)}

/* pills */
.pill{display:inline-block;padding:1px 8px;border-radius:9px;font-size:11px;
  letter-spacing:.04em;border:1px solid currentColor;font-weight:600}
.pill.gas{color:var(--ok)} .pill.marg{color:var(--warn)}
.pill.no{color:var(--bad)}

/* plates */
.plate{display:grid;grid-template-columns:minmax(0,1.05fr) minmax(0,1fr);
  gap:30px;align-items:start;padding:38px 0;border-bottom:1px solid var(--rule2)}
.plate:last-of-type{border-bottom:0}
@media(max-width:900px){.plate{grid-template-columns:1fr}}
.ptitle{display:flex;align-items:baseline;gap:12px;flex-wrap:wrap;
  margin:0 0 3px}
.ptitle h3{font-size:29px;margin:0;font-weight:600;letter-spacing:-.015em}
.ptitle .pr{color:var(--dim);font-size:14.5px}
.pmeta{color:var(--faint);font-size:12.5px;margin:0 0 16px}

/* viewer */
.viewer{border:1px solid var(--rule);border-radius:3px;background:var(--panel);
  overflow:hidden;box-shadow:var(--shadow)}
.vcanvas{display:block;width:100%;height:392px;cursor:grab;touch-action:none}
.vcanvas:active{cursor:grabbing}
.vbar{display:flex;flex-wrap:wrap;gap:6px;padding:9px 10px;
  border-top:1px solid var(--rule);align-items:center}
.vbar .sp{flex:1}
button.t{font:inherit;font-size:11.5px;letter-spacing:.03em;padding:3px 9px;
  border:1px solid var(--rule);background:transparent;color:var(--dim);
  border-radius:2px;cursor:pointer}
button.t:hover{border-color:var(--faint);color:var(--ink)}
button.t[aria-pressed="true"]{background:var(--six);border-color:var(--six);
  color:#fff}
button.t[aria-pressed="true"]{color:var(--onsix)}
.vlegend{display:flex;flex-wrap:wrap;gap:12px;padding:0 10px 10px;
  font-size:11.5px;color:var(--faint);align-items:center}
.sw{display:inline-block;width:11px;height:11px;border-radius:50%;
  vertical-align:-1px;margin-right:4px}
.swl{display:inline-block;width:15px;height:3px;vertical-align:2px;
  margin-right:4px;border-radius:2px}

/* data blocks */
.blocks{display:flex;flex-direction:column;gap:20px}
.block h4{font-size:10.5px;letter-spacing:.14em;text-transform:uppercase;
  color:var(--faint);margin:0 0 8px;font-weight:600;
  padding-bottom:5px;border-bottom:1px solid var(--rule)}
.kv{display:grid;grid-template-columns:auto 1fr;gap:3px 16px;font-size:13.5px}
.kv dt{color:var(--dim)}
.kv dd{margin:0;text-align:right;font-family:ui-monospace,SFMono-Regular,Menlo,
  monospace;font-variant-numeric:tabular-nums}
.orb{width:100%;font-size:12.5px}
.orb td,.orb th{padding:4px 7px}
.bar{display:flex;height:9px;border-radius:2px;overflow:hidden;margin-top:7px;
  border:1px solid var(--rule)}
.bar i{display:block}
.zkey{display:flex;gap:11px;flex-wrap:wrap;font-size:11.5px;
  color:var(--faint);margin-top:6px}
.callout{border-left:3px solid var(--six);padding:2px 0 2px 13px;
  font-size:13.5px;color:var(--dim);margin-top:4px}
.callout b{color:var(--ink);font-weight:600}
footer{padding:34px 0 60px;color:var(--faint);font-size:12.5px}
footer p{max-width:74ch}
.tag{display:inline-block;font-size:10.5px;letter-spacing:.1em;
  text-transform:uppercase;color:var(--warn);border:1px solid currentColor;
  padding:1px 7px;border-radius:2px}
"""

JS = r"""
(function(){
'use strict';
var DATA = window.__LIB__, RM = matchMedia('(prefers-reduced-motion: reduce)');
function css(n){return getComputedStyle(document.documentElement)
  .getPropertyValue(n).trim();}

function Viewer(host, sets){
  var cv = host.querySelector('canvas'), ctx = cv.getContext('2d');
  var st = {set:'chunk', mode:'all', atoms:true, spin:!RM.matches,
            rx:-0.42, ry:0.62, zoom:1, vis:false, dirty:true};
  var M = new Float32Array(9);

  function rot(){
    var cx=Math.cos(st.rx), sx=Math.sin(st.rx),
        cy=Math.cos(st.ry), sy=Math.sin(st.ry);
    // Ry then Rx
    M[0]=cy;    M[1]=0;  M[2]=-sy;
    M[3]=sx*sy; M[4]=cx; M[5]=sx*cy;
    M[6]=cx*sy; M[7]=-sx;M[8]=cx*cy;
  }

  function draw(){
    var S = sets[st.set];
    var dpr = Math.min(2, window.devicePixelRatio||1);
    var W = cv.clientWidth, H = cv.clientHeight;
    if(cv.width !== (W*dpr|0) || cv.height !== (H*dpr|0)){
      cv.width = W*dpr|0; cv.height = H*dpr|0;
    }
    ctx.setTransform(dpr,0,0,dpr,0,0);
    ctx.clearRect(0,0,W,H);
    rot();
    var n=S.n, p=S.p, X=S.X, Y=S.Y, Z=S.Z;
    var s = Math.min(W,H)*0.49*st.zoom/1000;
    var ox=W/2, oy=H/2, zmin=1e9, zmax=-1e9;
    for(var i=0;i<n;i++){
      var a=p[3*i], b=p[3*i+1], c=p[3*i+2];
      var x=M[0]*a+M[1]*b+M[2]*c, y=M[3]*a+M[4]*b+M[5]*c,
          z=M[6]*a+M[7]*b+M[8]*c;
      X[i]=ox+x*s; Y[i]=oy-y*s; Z[i]=z;
      if(z<zmin)zmin=z; if(z>zmax)zmax=z;
    }
    var span=(zmax-zmin)||1, NB=5;
    var six=css('--six'), five=css('--five');
    if(S.bp){                                    // unit-cell frame
      for(var q0=0;q0<S.bn;q0++){
        var ba=S.bp[3*q0], bb=S.bp[3*q0+1], bc=S.bp[3*q0+2];
        S.bX[q0]=ox+(M[0]*ba+M[1]*bb+M[2]*bc)*s;
        S.bY[q0]=oy-(M[3]*ba+M[4]*bb+M[5]*bc)*s;
      }
      ctx.beginPath();
      for(var q1=0;q1<S.be.length;q1+=2){
        ctx.moveTo(S.bX[S.be[q1]],S.bY[S.be[q1]]);
        ctx.lineTo(S.bX[S.be[q1+1]],S.bY[S.be[q1+1]]);
      }
      ctx.globalAlpha=0.5; ctx.lineWidth=1;
      ctx.setLineDash([4,4]); ctx.strokeStyle=css('--faint');
      ctx.stroke(); ctx.setLineDash([]);
    }
    // edges, bucketed by depth so far bonds recede
    var e=S.e, d=S.d, ne=d.length;
    for(var cls=0; cls<2; cls++){
      var want = cls===0?5:6;
      if(st.mode==='six' && want===5) continue;
      if(st.mode==='five'&& want===6) continue;
      for(var band=0; band<NB; band++){
        ctx.beginPath();
        var any=false;
        for(var k=0;k<ne;k++){
          if((d[k]>=6?6:5)!==want) continue;
          var u=e[2*k], v=e[2*k+1];
          var zm=(Z[u]+Z[v])*0.5;
          var bb=((zm-zmin)/span*NB)|0; if(bb>=NB)bb=NB-1;
          if(bb!==band) continue;
          ctx.moveTo(X[u],Y[u]); ctx.lineTo(X[v],Y[v]); any=true;
        }
        if(!any) continue;
        var t=(band+0.5)/NB;                       // 0 = far, 1 = near
        ctx.globalAlpha = want===6 ? (0.30+0.70*t) : (0.12+0.42*t);
        ctx.lineWidth   = want===6 ? (0.9+1.5*t)   : (0.5+0.5*t);
        ctx.strokeStyle = want===6 ? six : five;
        ctx.stroke();
      }
    }
    // atoms
    if(st.atoms){
      var zc=[css('--z12'),css('--z14'),css('--z15'),css('--z16')];
      var zz=S.z, core=S.ncore||0;
      var idx=S.order;
      for(var q=0;q<n;q++){
        var i2=idx[q];
        var t2=(Z[i2]-zmin)/span;
        var r=(1.7+2.2*t2)*Math.sqrt(st.zoom);
        var kz=zz[i2], ci=kz>=16?3:kz>=15?2:kz>=14?1:0;
        ctx.globalAlpha=0.35+0.65*t2;
        ctx.fillStyle=zc[ci];
        ctx.beginPath(); ctx.arc(X[i2],Y[i2],r,0,6.2832); ctx.fill();
        if(core && i2<core){
          ctx.globalAlpha=0.9; ctx.lineWidth=1;
          ctx.strokeStyle=css('--six'); ctx.stroke();
        }
      }
    }
    ctx.globalAlpha=1;
  }

  function tick(){
    if(st.vis && st.spin && !RM.matches){ st.ry += 0.0032; st.dirty=true; }
    if(st.dirty && st.vis){ st.dirty=false; draw(); }
    requestAnimationFrame(tick);
  }

  // interaction
  var drag=null;
  cv.addEventListener('pointerdown',function(ev){
    drag={x:ev.clientX,y:ev.clientY}; st.spin=false;
    host.querySelector('[data-a="spin"]').setAttribute('aria-pressed','false');
    cv.setPointerCapture(ev.pointerId);
  });
  cv.addEventListener('pointermove',function(ev){
    if(!drag) return;
    st.ry += (ev.clientX-drag.x)*0.0075;
    st.rx += (ev.clientY-drag.y)*0.0075;
    st.rx = Math.max(-1.5,Math.min(1.5,st.rx));
    drag={x:ev.clientX,y:ev.clientY}; st.dirty=true;
  });
  cv.addEventListener('pointerup',function(){drag=null;});
  cv.addEventListener('pointercancel',function(){drag=null;});
  cv.addEventListener('wheel',function(ev){
    ev.preventDefault();
    st.zoom *= Math.exp(-ev.deltaY*0.0012);
    st.zoom = Math.max(0.35,Math.min(6,st.zoom));
    st.dirty=true;
  },{passive:false});

  host.addEventListener('click',function(ev){
    var b=ev.target.closest('button.t'); if(!b) return;
    var a=b.dataset.a, v=b.dataset.v;
    if(a==='set'||a==='mode'){
      st[a]=v;
      host.querySelectorAll('button[data-a="'+a+'"]').forEach(function(o){
        o.setAttribute('aria-pressed', String(o.dataset.v===v));});
    } else {
      st[a]=!st[a]; b.setAttribute('aria-pressed',String(st[a]));
    }
    st.dirty=true;
  });

  new IntersectionObserver(function(en){
    st.vis=en[0].isIntersecting; if(st.vis) st.dirty=true;
  },{rootMargin:'160px'}).observe(cv);
  addEventListener('resize',function(){st.dirty=true;});
  new MutationObserver(function(){st.dirty=true;})
    .observe(document.documentElement,{attributes:true,
      attributeFilter:['data-theme']});
  matchMedia('(prefers-color-scheme:dark)')
    .addEventListener('change',function(){st.dirty=true;});
  tick();
}

function prep(S){
  var n=S.n;
  S.p=Float32Array.from(S.p); S.e=Int32Array.from(S.e);
  S.d=Int8Array.from(S.d); S.z=Int16Array.from(S.z);
  S.X=new Float32Array(n); S.Y=new Float32Array(n); S.Z=new Float32Array(n);
  if(S.bp){ S.bp=Float32Array.from(S.bp); S.be=Int32Array.from(S.be);
            S.bn=S.bp.length/3;
            S.bX=new Float32Array(S.bn); S.bY=new Float32Array(S.bn); }
  // painter order for atoms: far to near, recomputed cheaply per frame is
  // wasteful, so sort once by index and rely on depth-scaled alpha instead
  S.order=new Int32Array(n); for(var i=0;i<n;i++) S.order[i]=i;
  return S;
}

document.querySelectorAll('[data-crystal]').forEach(function(host){
  var k=host.dataset.crystal, c=DATA.crystals[k];
  Viewer(host, {cell:prep(c.scene_cell), chunk:prep(c.scene_chunk)});
});
})();
"""


def fmt(x, n=3):
    return f"{x:.{n}f}"


def orbit_rows(orb):
    """(kind label, total, [(decoration, count)]) in a fixed display order."""
    kinds = [("vertex", "vertices", ZC),
             ("edge", "edges", ["deg5", "deg6"]),
             ("face", "triangles", ["555", "556", "566", "666"]),
             ("tet", "tets", None)]
    out = []
    for key, lbl, pref in kinds:
        by = orb[key]["by"]
        keys = ([k for k in pref if k in by] + sorted(k for k in by
                                                      if k not in (pref or []))
                ) if pref else sorted(by, key=lambda s: (s.count("6"), s))
        out.append((lbl, orb[key]["total"], [(k, by[k]) for k in keys]))
    return out


def build(lib, out_path):
    C = lib["crystals"]
    ef = lib["e_flat"]
    P = []
    A = P.append

    A('<title>The TCP crystal library</title>')
    A(f"<style>{CSS}</style>")

    # ---------------- masthead
    A('<header class="mast"><div class="wrap">')
    A('<p class="eyebrow">Discrete differential geometry · reference crystals</p>')
    A('<h1 class="serif">Ten tetrahedrally close-packed crystals, '
      'measured against flatness</h1>')
    A('<p class="lede">Every structure here is a triangulation of the '
      '3-torus in which each edge is shared by exactly five or six tetrahedra. '
      'That constraint fixes almost everything else — and it leaves one number '
      'free per crystal: the mean edge degree. This catalog sets all ten '
      'against the value that would make the Regge curvature vanish on '
      'average.</p>')
    A('<dl class="thesis">')
    A('<div><dt>Flat edge degree</dt><dd class="mono serif">e* = 2π ⁄ arccos(⅓)'
      '</dd><div class="sm mono">= ' + f"{ef:.7f}" + '</div></div>')
    A('<div><dt>Crystals</dt><dd class="mono">10</dd>'
      '<div class="sm">a15 · σ · Z · P · δ · R · μ · C15 · C14 · C36</div></div>')
    A('<div><dt>Closest to flat</dt><dd class="mono">R</dd>'
      '<div class="sm mono">e − e* = ' + f"{C['r']['e_native']-ef:+.7f}" +
      '</div></div>')
    A('<div><dt>Furthest</dt><dd class="mono">A15</dd>'
      '<div class="sm mono">e − e* = ' + f"{C['a15']['e_native']-ef:+.7f}" +
      '</div></div>')
    A('</dl></div></header>')

    # ---------------- overview
    A('<section><div class="wrap">')
    A('<h2>Index</h2><p class="h2sub serif">The library at a glance</p>')
    A('<p class="note">Ordered by mean edge degree, from the most curved-'
      'positive down through flat. The sign of e − e* is not cosmetic: it '
      'decides which Pachner move a crystal must use to reach flatness, and '
      'therefore what its defects are made of.</p>')
    A('<div class="scroll"><table><thead><tr>'
      '<th>Crystal</th><th class="n">Atoms<br>/cell</th>'
      '<th class="n">⟨CN⟩</th><th class="n">e native</th>'
      '<th class="n">e − e*</th>'
      '<th class="n grp">Z12</th><th class="n">Z14</th>'
      '<th class="n">Z15</th><th class="n">Z16</th>'
      '<th class="n grp">|Aut|/m³</th><th class="n">orbits<br>0/1/2/3</th>'
      '<th class="n grp">gas</th><th class="n">defect %</th>'
      '</tr></thead><tbody>')
    for k in ORDER:
        c = C[k]
        de = c["e_native"] - ef
        g = c["gas"]
        o = c["orb"]
        if g is None:
            verd, vtxt, dfrac = "no", "—", "—"
        elif k == "a15":
            verd, vtxt, dfrac = "marg", "marginal", f'{100*g["dfrac"]:.1f}'
        else:
            verd, vtxt, dfrac = "gas", "yes", f'{100*g["dfrac"]:.1f}'
        cls = ' class="sep"' if k == "r" else ""
        cells = [f'<td class="name">{NICE[k]}'
                 f'<span class="p">{html.escape(c["proto"])}</span></td>',
                 f'<td class="n">{c["ns"]}</td>',
                 f'<td class="n">{c["cn"]:.4f}</td>',
                 f'<td class="n">{c["e_native"]:.7f}</td>',
                 f'<td class="n">{de:+.6f}</td>']
        cells += [f'<td class="n{" grp" if j == 0 else ""}">'
                  f'{c["census"].get(z, 0) or "·"}</td>'
                  for j, z in enumerate(ZC)]
        cells += [f'<td class="n grp">{c["point"]}</td>',
                  f'<td class="n">{o["vertex"]["total"]}/{o["edge"]["total"]}/'
                  f'{o["face"]["total"]}/{o["tet"]["total"]}</td>',
                  f'<td class="n grp"><span class="pill {verd}">{vtxt}</span></td>',
                  f'<td class="n">{dfrac}</td>']
        A(f"<tr{cls}>" + "".join(cells) + "</tr>")
    A('</tbody></table></div>')
    A('<p class="note" style="margin-top:16px">R sits one part in 10⁵ from '
      'flat — the horizontal rule marks where e − e* changes sign. Above it, '
      'crystals are too curved and must <em>remove</em> edge share; below, too '
      'flat-poor and must add it.</p>')
    A('</div></section>')

    # ---------------- plates
    A('<section><div class="wrap">')
    A('<h2>Plates</h2><p class="h2sub serif">One crystal at a time</p>')
    A('<p class="note">Drag to rotate, scroll to zoom. <b>Degree-6 edges are '
      'drawn heavy and amber</b> — they are the disclination network, the '
      'skeleton along which the structure differs from an ideal tetrahedral '
      'packing. The unit-cell view outlines the cell\'s own atoms in amber and '
      'draws every one of their coordination shells whole.</p>')

    for k in ORDER:
        c = C[k]
        g = c["gas"]
        de = c["e_native"] - ef
        cell = c["cell"]
        hexish = abs(cell["gamma"] - 120) < 1e-6
        A(f'<div class="plate" data-crystal="{k}">')
        # --- left: viewer
        A('<div>')
        A(f'<div class="ptitle"><h3 class="serif">{NICE[k]}</h3>'
          f'<span class="pr">{html.escape(c["proto"])}</span></div>')
        A(f'<p class="pmeta mono">{c["ns"]} atoms/cell · block '
          f'{c["K"]}×{c["K"]}×{c["K"]} = {c["scene_chunk"]["n"]} atoms · '
          f'{len(c["scene_chunk"]["d"])} bonds drawn</p>')
        A('<div class="viewer">')
        A('<canvas class="vcanvas" role="img" aria-label="Rotatable 3D model '
          f'of the {html.escape(NICE[k])} structure with degree-6 edges '
          'highlighted"></canvas>')
        A('<div class="vbar">'
          '<button class="t" data-a="set" data-v="chunk" aria-pressed="true">'
          'Block</button>'
          '<button class="t" data-a="set" data-v="cell" aria-pressed="false">'
          'Unit cell</button><span class="sp"></span>'
          '<button class="t" data-a="mode" data-v="all" aria-pressed="true">'
          'All bonds</button>'
          '<button class="t" data-a="mode" data-v="six" aria-pressed="false">'
          'Degree 6 only</button>'
          '<button class="t" data-a="mode" data-v="five" aria-pressed="false">'
          'Degree 5 only</button>'
          '<button class="t" data-a="atoms" aria-pressed="true">Atoms</button>'
          '<button class="t" data-a="spin" aria-pressed="true">Spin</button>'
          '</div>')
        A('<div class="vlegend">'
          '<span><i class="swl" style="background:var(--six)"></i>degree 6 '
          '(disclination)</span>'
          '<span><i class="swl" style="background:var(--five)"></i>degree 5'
          '</span>'
          + "".join(f'<span><i class="sw" style="background:var(--z{z[1:]})">'
                    f'</i>{z}</span>' for z in ZC if c["census"].get(z))
          + '</div>')
        A('</div></div>')

        # --- right: data
        A('<div class="blocks">')
        A('<div class="block"><h4>Cell &amp; curvature</h4><dl class="kv">')
        lat = (f'a={cell["a"]:.4f}  c={cell["c"]:.4f}' if hexish
               else f'{cell["a"]:.4f} × {cell["b"]:.4f} × {cell["c"]:.4f}')
        A(f'<dt>Lattice (a=1 units)</dt><dd>{lat}</dd>')
        A(f'<dt>Angles α β γ</dt><dd>{cell["alpha"]:.1f}° {cell["beta"]:.1f}° '
          f'{cell["gamma"]:.1f}°</dd>')
        A(f'<dt>Atoms per cell</dt><dd>{c["ns"]}</dd>')
        A(f'<dt>Mean coordination ⟨CN⟩</dt><dd>{c["cn"]:.5f}</dd>')
        A(f'<dt>Mean edge degree ⟨e⟩</dt><dd>{c["e_native"]:.7f}</dd>')
        A(f'<dt>e − e*</dt><dd>{de:+.7f}</dd>')
        A(f'<dt>Pin gap at m={c["mcell"]}</dt><dd>{c["gap"]:+.2f} f₁</dd>')
        A(f'<dt>Forced defect debt</dt><dd>{c["nforced"]:.0f} moves '
          f'({1000*c["nforced"]/c["fv"][0]:.0f}/1000 v)</dd>')
        A('</dl>')
        tot = sum(c["census"].values())
        A('<div class="bar">' + "".join(
            f'<i style="flex:{c["census"][z]};background:var(--z{z[1:]})"></i>'
            for z in ZC if c["census"].get(z)) + '</div>')
        A('<div class="zkey">' + " ".join(
            f'<span>{z} × {c["census"][z]}</span>'
            for z in ZC if c["census"].get(z))
          + f'<span>({tot} atoms)</span></div>')
        A('</div>')

        A(f'<div class="block"><h4>Symmetry — exact Aut(K)</h4>'
          f'<dl class="kv">'
          f'<dt>|Aut| at m={c["mcell"]}</dt><dd>{c["aut"]}</dd>'
          f'<dt>|Aut| ⁄ m³ (space group order)</dt><dd>{c["point"]}</dd>'
          f'</dl>')
        A('<table class="orb"><thead><tr><th>Simplices</th>'
          '<th class="n">Orbits</th><th>Refined by decoration</th>'
          '</tr></thead><tbody>')
        for lbl, total, items in orbit_rows(c["orb"]):
            det = " · ".join(f'{html.escape(d)} <b>{n}</b>' for d, n in items)
            A(f'<tr><td>{lbl}</td><td class="n">{total}</td>'
              f'<td style="text-align:left;white-space:normal" class="mono">'
              f'{det}</td></tr>')
        A('</tbody></table></div>')

        if g:
            note = {"paid": "discharges the pin debt",
                    "unpaid": "does <b>not</b> discharge the pin debt",
                    "n/a": "no debt to discharge — already flat"}[g["pin_note"]]
            A('<div class="block"><h4>Defect gas at c<sub>imp</sub> = '
              f'{g["cimp"]:.2f}</h4><dl class="kv">'
              f'<dt>Defect vertices</dt><dd>{g["n_ill"]:.0f} '
              f'({100*g["dfrac"]:.1f}%)</dd>'
              f'<dt>Complexes</dt><dd>{g["ncomp"]:.0f}</dd>'
              f'<dt>Largest complex</dt><dd>{g["top1"]:.1f} vertices</dd>'
              f'<dt>Max impurity valence</dt><dd>{g["max_m"]}</dd>'
              f'<dt>Turnover / 50 sweeps</dt><dd>{g["turnover"]:.2f}</dd>'
              f'<dt>Drift significance</dt><dd>{g["drift_z"]:.1f} σ</dd>'
              f'<dt>Achieved ⟨e⟩</dt><dd>{g["e_mean"]:.6f}</dd>'
              '</dl>'
              f'<p class="callout">The gas {note}.'
              + (' <span class="tag">not equilibrated</span> This host\'s '
                 'population was still climbing at 6000 sweeps (15 → 66); '
                 'treat the verdict as marginal.' if k == "a15" else "")
              + '</p></div>')
        A('</div></div>')
    A('</div></section>')

    # ---------------- findings
    A('<section><div class="wrap">')
    A('<h2>Across the library</h2>'
      '<p class="h2sub serif">Three things the orbit tables settle</p>')
    A('<div style="display:grid;grid-template-columns:repeat(auto-fit,'
      'minmax(270px,1fr));gap:26px;margin-top:8px">')
    A('<div><h4 style="font-size:14px;margin:0 0 7px">No triangle carries two '
      'degree-6 edges</h4><p class="note">Only (5,5,5) and (5,5,6) triples '
      'occur — in all ten crystals, every orbit. Since any two edges of a '
      'tetrahedron that share a vertex also span a face, this forces the tet '
      'table: a cell\'s degree-6 subgraph is empty, a single edge, or a pair '
      'of <em>opposite</em> edges. Nothing denser exists. Switch a viewer to '
      '“Degree 6 only” and the network has no corners.</p></div>')
    A('<div><h4 style="font-size:14px;margin:0 0 7px">Aut(K) is exactly the '
      'space group</h4><p class="note">|Aut| ⁄ m³ reproduces the '
      'crystallographic order for every structure, and vertex orbits match '
      'Wyckoff site counts. No accidental combinatorial symmetry — worth '
      'checking, because the automorphism group of an abstract complex is '
      'free to exceed the physical one. Counts are supercell-independent: '
      'recomputing at m = 3 changes |Aut| but not one orbit count.</p></div>')
    A('<div><h4 style="font-size:14px;margin:0 0 7px">δ is the only free '
      'action</h4><p class="note">Every δ stabilizer is trivial, so its orbit '
      'counts are just f<sub>k</sub> ⁄ |Aut| — 14, 94, 160, 80 — and its '
      'alternating sum vanishes automatically. The other nine have orbits with '
      'nontrivial stabilizers, so only the stabilizer-weighted sum is pinned '
      'to χ ⁄ |Aut| = 0.</p></div>')
    A('</div>')
    A('<p class="note" style="margin-top:30px"><b>P and δ share an f-vector '
      'exactly</b> (1512, 10152, 17280, 8640) and the same mean coordination, '
      'yet differ by a factor two in symmetry order — 216 against 108 — and by '
      '66 tet orbits. Identical counting, different crystals.</p>')
    A('</div></section>')

    A('<footer><div class="wrap">')
    A('<p>Structures built from published Wyckoff positions by '
      '<span class="mono">scripts/tcp_reference.py</span> and validated on '
      'construction: expected f-vector, no duplicate facets, Euler '
      'characteristic 0, orientable, the exact link sum rule, all edge degrees '
      'in {5,6}, and the literature Z-class census. Symmetry from '
      '<span class="mono">symmetry.CrystalSymmetry</span> (exact automorphisms '
      'as explicit permutations, not Weisfeiler–Leman colours). Gas figures '
      'are single-chain, 1000-sweep burn plus 2000 sampled, stationarity by '
      'late-window slope against a moving-block bootstrap σ — '
      '<b>provisional: no two-sided R-hat certification</b>.</p>')
    A('<p style="margin-top:12px">Geometry drawn here is the periodic Delaunay '
      'triangulation quotiented to T³; hinge degrees are computed on the '
      'quotient and are therefore exact bulk values. In the block view, bonds '
      'that wrap the torus are omitted, so the outer surface is cut.</p>')
    A('</div></footer>')

    A('<script>window.__LIB__=' + json.dumps(lib, separators=(",", ":"))
      + ';</script>')
    A(f'<script>{JS}</script>')

    with open(out_path, "w") as fh:
        fh.write("\n".join(P))
    return out_path


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", default=os.path.join(
        _ROOT, "data", "crystal_gas", "library.json"))
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    lib = json.load(open(args.data))
    p = build(lib, args.out)
    print(f"-> {p}  ({os.path.getsize(p)/1e6:.2f} MB)")


if __name__ == "__main__":
    main()
