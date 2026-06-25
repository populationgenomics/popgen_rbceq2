# Writing your repo's PRODUCT.md (north star)

Every durable repo needs a north-star doc that an agent and a new contributor
read first — what the thing is, why it exists, and the load-bearing decisions
that shape it. This file tells you how to write your own; **replace it with the
doc you write.**

## What this doc is (and isn't)

- It's the **frame to start detailed plans from**, so context isn't re-litigated
  every session. Keep it current, but this is not expected to change frequently.
- It is **not the build plan** — concrete current work lives in issues / design
  docs / plans; keep those out of here.
- Write it **for an LLM reader**: dense, each fact stated **once**, no rhetorical
  emphasis or recap. Every token is re-paid on every future read.

## What to capture

Cover these, adapting the headings to your project:

1. **What it is.** One tight statement of the thing, then its inputs, the
   structure it works within, and what it outputs.
2. **Why it exists.** The problem / bottleneck it attacks, why now, and what
   success looks like in observable terms.
3. **Who it's for.** The users — and whether they're also the people building it.
4. **The core thesis.** The actual durable contribution — the bet about where the
   lasting value sits, versus what's thin and replaceable.
5. **Load-bearing design principles.** The non-obvious choices everything else
   must respect (the ones you'd otherwise re-argue in every design review).
6. **Scope boundaries & ecosystem.** What it is **not**; the neighbouring
   tools/systems and exactly how this relates to each.
7. **Relevant external repositories and resources.** Point at the external repos
   an agent should *read* for this work — upstream tools, libraries, reference
   implementations, related pipelines — with a one-line note on what each is for,
   so it doesn't guess at code it could read. And where relevant context should
   be extracted or summarised from an external resource (e.g. documentation for a
   tool like DRAGEN or Hail), capture that summary here.
8. **Bets & open questions.** The bets you hold deliberately — each with *what
   would falsify it* — and the questions you'll resolve by building, not up front.
9. **The current slice.** The one narrowed target you're building first, with a
   pointer to where its concrete definition lives (not inlined here).

## Domain terms

Define the terms a newcomer (human or agent) would otherwise guess at in a
companion `GLOSSARY.md` (see its template), and point at it from here. Don't
inline a long term list in this doc.
