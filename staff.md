---
layout: page
title: Instructors
permalink: /instructors/
---
<!----- <h1>Instructors</h1> ---->


  {% for author in site.authors %}
  <section>
    <img itemprop="image" class="img-rounded" src="{{ site.baseurl }}/assets/img/{{ author.picture }}" alt="{{ author.name }}" height=150 width=150>
    <h2>{{ author.name }}</h2>
    <h4>{{ author.organization }}</h4>
    <h4>{{ author.position }}</h4>
    <p>{{ author.content | markdownify }}</p>
    <line></line>
  <br><br>
  {% endfor %}