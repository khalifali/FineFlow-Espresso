import fitz,json,zipfile
from pathlib import Path
root=Path(__file__).resolve().parent;deck=fitz.open(root/'fineflow_espresso_beamer.pdf');notes=json.loads((root/'transcript.json').read_text());out=fitz.open()
for i,s in enumerate(deck):
 p=out.new_page(width=1200,height=675)
 p.insert_text((28,40),'FineFlow-Espresso | Presenter transcript',fontsize=20,color=(0,.396,.741))
 p.insert_text((28,65),f'Slide {i+1} of {len(deck)}',fontsize=12,color=(.4,.4,.4))
 p.show_pdf_page(fitz.Rect(20,110,720,504),deck,i)
 p.draw_line((740,96),(740,615),color=(.8,.8,.8),width=.6)
 p.insert_text((770,118),'Presenter transcript',fontsize=17,color=(0,.396,.741))
 text=notes[str(i+1)]
 spare=p.insert_textbox(fitz.Rect(770,146,1170,615),text,fontsize=16,lineheight=1.45,fontname='helv',color=(.1,.1,.1))
 assert spare>=0,(i,spare)
 p.insert_text((28,647),'Revised discussion of Heiko Briesen\'s comments | 11 September 2026',fontsize=10,color=(.4,.4,.4))
out.save(root/'fineflow_espresso_with_transcript.pdf',garbage=4,deflate=True)
(root/'presenter_transcript.txt').write_text('\n\n'.join(f'SLIDE {i}\n{notes[str(i)]}' for i in range(1,len(deck)+1)))
with zipfile.ZipFile(root/'fineflow_espresso_source.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in [root/'fineflow_espresso_beamer.tex',root/'presenter_transcript.txt',*sorted((root/'figures').iterdir())]:z.write(p,p.relative_to(root))

print('Slides and transcript pages rebuilt; source archive ready')
