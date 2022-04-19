---
permalink: /instructors/
layout: default
title: Instructors
---
<h1>Instructors</h1>


  {% for author in site.authors %}
  <section>
    <img itemprop="image" class="img-rounded" src="{{ site.baseurl }}/assets/img/{{ author.picture }}" alt="{{ author.name }}" height=150 width=150>
    <h2>{{ author.name }}</h2>
    <h3>{{ author.organization }}</h3>
    <h3>{{ author.position }}</h3>
    <p>{{ author.content | markdownify }}</p>
    <line></line>
  <br>
  {% endfor %}
